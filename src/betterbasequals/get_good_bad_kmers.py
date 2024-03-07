import pysam
import py2bit
import os
from collections import Counter
from betterbasequals.utils import reverse_complement, Read, zip_pileups_single_chrom, eprint, open_bam_w_index, SortedBed, ostrand

# def get_mut_type(ref, alt):
#     if ref in ['T', 'G']:
#         mtype = reverse_complement(ref) + '->' + reverse_complement(alt)
#     else:
#         mtype = ref + '->' + alt
#     return mtype

def get_mut_type(ref, alt, kmer, is_reverse):
    if is_reverse:
        mtype = reverse_complement(ref) + '->' + reverse_complement(alt)
    else:
        mtype = ref + '->' + alt
    return mtype


def get_BQ_pair(BQ1, BQ2):
    if BQ1 > BQ2:
        return f'({BQ1},{BQ2})'
    else:
        return f'({BQ2},{BQ1})'

class MutationCounterWFilter:
    def __init__(
        self,
        bam_file,
        twobit_file,
        filter_bam_file,
        mapq = 50,
        min_base_qual = 1,
        filter_mapq = 20,
        min_base_qual_filter=20, 
        min_depth=1, 
        max_depth=1000000, 
        radius=3, 
        prefix="", 
        max_alt_frac = None,
        min_filter_count = 2,
        min_filter_depth = 10,
        max_filter_depth = 100,
        min_enddist = 6,
        max_mismatch = 2,
        bed_file = None,
        #output,
    ):
        self.bam_file = open_bam_w_index(bam_file)
        if filter_bam_file is None:
            self.filter_bam_file = None
            if max_alt_frac is None:
                #Set default to filter heterozygous variants if no filter is used:
                self.max_alt_frac = 0.2
        else:
            self.filter_bam_file = open_bam_w_index(filter_bam_file)
            if max_alt_frac is None:
                #Should I have different default here?
                #Maybe 1.0 corresponding to no filter?
                self.max_alt_frac = 0.2

        self.tb = py2bit.open(twobit_file)
   
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix
        self.min_filter_count = min_filter_count
        self.min_filter_depth = min_filter_depth
        self.max_filter_depth = max_filter_depth
        self.min_enddist = min_enddist
        self.max_mismatch = max_mismatch
        if not bed_file is None:
            self.bed = SortedBed(bed_file)
        else:
            self.bed = None

    def __del__(self):
        self.tb.close()

    def count_mutations_all_chroms(self):
        event_kmers = Counter()
        for idx_stats in self.bam_file.get_index_statistics():
            if idx_stats.mapped > 0:
                event_c = self.count_mutations(idx_stats.contig, None, None)
                event_kmers += event_c
        return event_kmers


    def count_mutations(self, chrom, start, stop):
        event_kmers = Counter()

        pileup = self.bam_file.pileup(
            contig = chrom,
            start = start,
            stop = stop,
            truncate = True,
            max_depth = 1000000,
            min_mapping_quality = self.mapq,
            ignore_overlaps = False,
            flag_require = 0,  # No requirement
            flag_filter = 3848,
            min_base_quality = self.min_base_qual,
        )
        if not self.filter_bam_file is None:
            filter_pileup = self.filter_bam_file.pileup(
                contig = chrom,
                start = start,
                stop = stop,
                truncate = True,
                min_mapping_quality = self.filter_mapq,
                ignore_overlaps = False,
                flag_require = 2,  # proper paired
                flag_filter = 3848,
                min_base_quality = self.min_base_qual_filter,
            )
            for pileupcolumn, filter_pc in zip_pileups_single_chrom(pileup, filter_pileup):
                N_filter = 0
                n_alleles = Counter()
                for pread in filter_pc.pileups:
                    pos = pread.query_position
                    if not pos is None:
                        n_alleles[pread.alignment.query_sequence[pos]] += 1
                        N_filter += 1
                if N_filter < self.min_filter_depth or N_filter > self.max_filter_depth:
                    continue                

                # We can set filter so that we only considder sites where
                # Any non-ref allele is seen min_filter_count times or more
                if len(n_alleles) > 1:
                    major, second = n_alleles.most_common(2)
                    if second[1] >= self.min_filter_count:
                        continue
                else:
                    major = n_alleles.most_common(1)[0]

                self.handle_pileup(pileupcolumn, event_kmers, major[0])
        else:
            for pileupcolumn in pileup:
                self.handle_pileup(pileupcolumn, event_kmers)
        
        return event_kmers

    def handle_pileup(self, pileupcolumn, event_kmers, major_allele = None):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%10000 ==0:
            eprint(f"{chrom}:{ref_pos}")    

        if not self.bed is None and self.bed.query(chrom, ref_pos):
            return
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue
        if ref_pos - self.radius < 0:
            return 

        kmer = self.tb.sequence(self.prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
        
        if 'N' in kmer:
            return
        
        if len(kmer) != 2*self.radius + 1:
            return 0
        
        ref = kmer[self.radius]

        #If the major in the filter bam file does not match the ref we ignore the site.
        if major_allele is not None and major_allele != ref:
            return

        if ref not in "ATGC":
            return
        
        #if ref not in ['A', 'C']:
        kmer_other_stand = reverse_complement(kmer)

        reads_mem = {}
        event_list = []
        coverage = 0

        has_good = set()
        n_alleles = Counter()

        for pileup_read in pileupcolumn.pileups:
            # test for deletion at pileup
            #if pileup_read.is_del or pileup_read.is_refskip:
            #    continue
            coverage += 1
            #TODO: should consider what the right solution is if there is deletion at overlap

            # fetch read information
            read = Read(pileup_read)
            
            # # test if read is okay
            # if (
            #     read.allel not in "ATGC"
            #     or read.start is None
            #     or read.end is None
            #     #or read.NH != 1
            # ):
            #     continue

            # Look for read partner
            if read.query_name in reads_mem:
                # found partner process read pair
                mem_read = reads_mem.pop(read.query_name)

                if (not read.is_good(self.min_enddist, self.max_mismatch, self.mapq)) or (not mem_read.is_good(self.min_enddist, self.max_mismatch, self.mapq)):
                    continue

                if read.allel == mem_read.allel:

                    event_list.append(('good', read.allel, read.base_qual, read.is_reverse, read.isR1))
                    event_list.append(('good', mem_read.allel, mem_read.base_qual, mem_read.is_reverse, mem_read.isR1))

                    #Count good tuples two times: once for each strand!
                    #TODO: Why?
                    event_list.append(('good_tuple', read.allel, get_BQ_pair(read.base_qual, mem_read.base_qual), read.is_reverse, read.isR1))
                    event_list.append(('good_tuple', read.allel, get_BQ_pair(read.base_qual, mem_read.base_qual), mem_read.is_reverse, mem_read.isR1))

                    if max(read.base_qual, mem_read.base_qual) > 30:
                        has_good.add(read.allel)
                        n_alleles[read.allel] += 1

                else:
                    if read.allel != ref and mem_read.allel == ref:# and mem_read.base_qual > 30:
                        event_list.append(('bad', read.allel, read.base_qual, read.is_reverse, read.isR1))
                        event_list.append(('good', mem_read.allel, mem_read.base_qual, mem_read.is_reverse, mem_read.isR1))
                        event_list.append(('bad_tuple', read.allel, get_BQ_pair(read.base_qual, mem_read.base_qual), read.is_reverse, read.isR1))

                    elif mem_read.allel != ref and read.allel == ref:# and read.base_qual > 30:
                        event_list.append(('bad', mem_read.allel, mem_read.base_qual, mem_read.is_reverse, mem_read.isR1))
                        event_list.append(('good', read.allel, read.base_qual, read.is_reverse, read.isR1))
                        event_list.append(('bad_tuple', mem_read.allel, get_BQ_pair(read.base_qual, mem_read.base_qual), mem_read.is_reverse, mem_read.isR1))

            else:            
                reads_mem[read.query_name] = read
        
        # Handle reads without partner (ie. no overlap)
        for read in reads_mem.values():
            if read.base_qual > 30:
                n_alleles[read.allel] += 1
            if read.is_good(self.min_enddist, self.max_mismatch, self.mapq):
                event_list.append(('singleton', read.allel, read.base_qual, read.is_reverse, read.isR1))

        if coverage < self.min_depth or coverage > self.max_depth:
            return
        
        N = sum(n_alleles.values())

        # We ignore sites where matching reference allele aren't seen in an overlap
        # OBS: Have commented out. 2023.06.02
        # In files with litlle overap we maybe should not require this.
        #if ref not in has_good:
        #    return
                
        if len(n_alleles) > 1:
            major, second = n_alleles.most_common(2)
            # If we see a high fraction of an alternative allele we ignore the site.
            if second[1]/N >= self.max_alt_frac:
                return
        elif len(n_alleles) == 1:
            major = n_alleles.most_common(1)[0]
        else:
            #We require that we see a good quality allele at the site
            return

        #If the major allele is not the ref allele we ignore the site.
        if major[0] != ref:
            return
        

        for event_type, allele, base_qual, is_reverse, is_R1 in event_list:
            if is_reverse:
                mut_type = f'{ostrand[ref]}->{ostrand[allele]}'
                this_kmer = kmer_other_stand
            else:
                mut_type = f'{ref}->{allele}'
                this_kmer = kmer
            #TODO: should this be for all categories?
            #if event_type in ['singleton', 'good_tuple', 'bad_tuple'] and n_alleles[allele]>1:
            #    continue
            if n_alleles[allele]>1:
                continue
            event_kmers[(event_type, base_qual, mut_type, int(is_R1), this_kmer)] += 1

