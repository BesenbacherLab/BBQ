import pysam
import py2bit
import os
from collections import Counter
from betterbasequals.utils import reverse_complement, Read, zip_pileups_single_chrom, eprint, open_bam_w_index

def get_mut_type(ref, alt):
    if ref in ['T', 'G']:
        mtype = reverse_complement(ref) + '->' + reverse_complement(alt)
    else:
        mtype = ref + '->' + alt
    return mtype


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
        max_bad_frac=0.05,
        min_filter_count = 2,
        #output,
    ):
        self.bam_file = open_bam_w_index(bam_file)
        if filter_bam_file is None:
            self.filter_bam_file = None
        else:
            self.filter_bam_file = open_bam_w_index(filter_bam_file)
        self.tb = py2bit.open(twobit_file)
   
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix
        self.max_bad_alt = max_bad_frac
        self.min_filter_count = min_filter_count

    def __del__(self):
        self.tb.close()

    def count_mutation_all_chroms(self):
        good_kmers = Counter()
        bad_kmers = Counter()
        for idx_stats in self.bam_file.get_index_statistics():
            if idx_stats.mapped > 0:
                good_c, bad_c = self.count_mutations(idx_stats.contig, None, None)
                good_kmers += good_c
                bad_kmers += bad_c
        return good_kmers, bad_kmers


    def count_mutations(self, chrom, start, stop):
        #if self.base_quals is None:
        #    good_kmers = {y:{x:Counter() for x in mtypes} for y in range(1,99)}
        #    bad_kmers = {y:{x:Counter() for x in mtypes} for y in range(1,99)}
        #else:
        #    good_kmers = {y:{x:Counter() for x in mtypes} for y in self.base_quals}
        #    bad_kmers = {y:{x:Counter() for x in mtypes} for y in self.base_quals}
        good_kmers = Counter()
        bad_kmers = Counter()

        pileup = self.bam_file.pileup(
            contig = chrom,
            start = start,
            stop = stop,
            truncate = True,
            max_depth = 1000000,
            min_mapping_quality = self.mapq,
            ignore_overlaps = False,
            flag_require = 2,  # proper paired
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
                        n_alleles.add(pread.alignment.query_sequence[pos])
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
                    major = n_alleles.most_common(1)

                self.handle_pileup(pileupcolumn, good_kmers, bad_kmers, major[0])
        else:
            for pileupcolumn in pileup:
                self.handle_pileup(pileupcolumn, good_kmers, bad_kmers)
        
        return good_kmers, bad_kmers

    def handle_pileup(self, pileupcolumn, good_kmers, bad_kmers, major_allele = None):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%10000 ==0:
            eprint(f"{chrom}:{ref_pos}")            
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue
        if ref_pos-self.radius < 0:
            return 

        kmer = self.tb.sequence(self.prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
        if 'N' in kmer:
            return

        ref = kmer[self.radius]

        if major_allele is not None and major_allele != ref:
            return

        if ref not in "ATGC":
            return
        if ref not in ['A', 'C']:
            kmer = reverse_complement(kmer)

        reads_mem = {}
        event_list = []
        coverage = 0

        has_good = set()
        n_allele = Counter()

        for pileup_read in pileupcolumn.pileups:
            # test for deletion at pileup
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            coverage += 1
            #TODO: should consider what the right solution is if there is deletion at overlap

            # fetch read information
            read = Read(pileup_read)

            # test if read is okay
            if (
                read.allel not in "ATGC"
                or read.start is None
                or read.end is None
                #or read.NH != 1
            ):
                continue

            if read.base_qual > 30:
                n_allele[read.allel] += 1

            # Look for read partner
            if read.query_name in reads_mem:
                # found partner process read pair
                mem_read = reads_mem.pop(read.query_name)

                if read.allel == mem_read.allel:
                    event_list.append(('good', read.allel, read.base_qual))
                    event_list.append(('good', read.allel, mem_read.base_qual))
                    
                    if max(read.base_qual, mem_read.base_qual) > 30:
                        has_good.add(read.allel)

                else:
                    #Are ignoring ref alleles pt... should adjust later

                    if read.allel != ref:
                        event_list.append(('bad', read.allel, read.base_qual))

                    elif mem_read.allel != ref:
                        event_list.append(('bad', mem_read.allel, mem_read.base_qual))
            else:            
                reads_mem[read.query_name] = read
        if coverage < self.min_depth or coverage > self.max_depth:
            return
        N = sum(n_allele.values())

        #We ignore sites where matching reference allele aren't seen in an overlap
        if ref not in has_good:
            return
        
        #
        #if len(n_alleles) > 1:
        #    major, second = n_alleles.most_common(2)
        #    if second[1] >= self.min_filter_count:
        #        continue
        #    else:
        #        major = n_alleles.most_common(1)
        
        #if (n_allele[allele]/N) > self.max_bad_frac:
            likely_variant = True
        #if major_allele is None:
        #    

        for event_type, allele, base_qual in event_list:
            mut_type = get_mut_type(ref, allele)     
            if event_type == 'good':
                if ref != allele or len(has_good) == 1:
                    good_kmers[(base_qual, mut_type, kmer)] += 1
                #good_kmers[base_qual][mut_type][kmer] += 1
            elif event_type == 'bad':
                #OBS: Have commented this out. 2023.02.02
                # Isince it would create a bias to have a filter on bad sites that
                # is not applied to good sites. 
                # We do not considder bad variants at sites where
                #if (allele not in has_good and 
                #    (n_allele[allele]/N) < self.max_bad_frac):
                # 
                if len(has_good) == 1:
                    bad_kmers[(base_qual, mut_type, kmer)] += 1
                    #bad_kmers[base_qual][mut_type][kmer] += 1
