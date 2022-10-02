from asyncio import events
import pysam
import py2bit
import os
from collections import Counter
from betterbasequals.utils import reverse_complement, mtypes, Read, zip_pileups, eprint, open_bam_w_index

def get_mut_type(ref, alt):
    if ref in ['T', 'G']:
        mtype = reverse_complement(ref) + '->' + reverse_complement(alt)
    else:
        mtype = ref + '->' + alt
    return mtype


class MutationFinderWFilter:
    def __init__(
        self,
        bam_file,
        twobit_file,
        filter_bam_file,
        basequals = [11,25,37]
        #output,
    ):
        self.bam_file = open_bam_w_index(bam_file)
        if filter_bam_file is None:
            self.filter_bam_file = None
        else:
            self.filter_bam_file = open_bam_w_index(filter_bam_file)
        self.tb = py2bit.open(twobit_file)
        self.base_quals = basequals

    def __del__(self):
        self.tb.close()

    def find_mutations(self, chrom, start, stop, mapq=50, min_base_qual=1, filter_mapq =20, min_base_qual_filter=20, min_depth=1, max_depth=5000, radius=3, prefix="",):
        good_kmers = {y:{x:Counter() for x in mtypes} for y in self.base_quals}
        bad_kmers = {y:{x:Counter() for x in mtypes} for y in self.base_quals}
        pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
            min_base_quality = min_base_qual,
        )
        if not self.filter_bam_file is None:
            filter_pileup = self.filter_bam_file.pileup(
                contig=chrom,
                start=start,
                stop=stop,
                truncate=True,
                min_mapping_quality=filter_mapq,
                ignore_overlaps=False,
                flag_require=2,  # proper paired
                flag_filter=3848,
                min_base_quality = min_base_qual_filter,
            )
            for pileupcolumn, filter_pc in zip_pileups(pileup, filter_pileup):
                N_filter = 0
                for pread in filter_pc.pileups:
                    N_filter += 1
                #TODO: replace hardcoded numbers with values relative to mean coverage
                if N_filter < 25 or N_filter > 55:
                    continue
                self.handle_pileup(pileupcolumn, good_kmers, bad_kmers, prefix, min_depth, max_depth, radius)
        else:
            for pileupcolumn in pileup:
                self.handle_pileup(pileupcolumn, good_kmers, bad_kmers, prefix, min_depth, max_depth, radius)
        
        return good_kmers, bad_kmers

    def handle_pileup(self, pileupcolumn, good_kmers, bad_kmers, prefix, min_depth, max_depth, radius):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%10000 ==0:
            eprint(f"{chrom}:{ref_pos}")            
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue
        kmer = self.tb.sequence(prefix + chrom, ref_pos- radius, ref_pos + radius + 1)
        if 'N' in kmer:
            return
        ref = kmer[radius]
        if ref not in "ATGC":
            return

        reads_mem = dict()
        event_list = []
        coverage = 0
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
            # Look for read partner
            if read.query_name in reads_mem:
                # found partner process read pair
                mem_read = reads_mem.pop(read.query_name)

                # We ignore alleles where pair doesn't match (could add else clause to handle)
                if read.allel == mem_read.allel:
                    #Are ignoring ref alleles pt... should adjust later
                    if read.allel == ref:
                        continue

                    mut_type = get_mut_type(ref, read.allel)     
                    event_list.append(('good', mut_type, read.base_qual))
                    event_list.append(('good', mut_type, mem_read.base_qual))

                else:
                    #Are ignoring ref alleles pt... should adjust later
                    if read.allel != ref:
                        mut_type = get_mut_type(ref, read.allel)
                        event_list.append(('bad', mut_type, read.base_qual))

                    elif mem_read.allel != ref:
                        mut_type = get_mut_type(ref, mem_read.allel)
                        event_list.append(('bad', mut_type, mem_read.base_qual))
            else:            
                reads_mem[read.query_name] = read
        if coverage < min_depth or coverage > max_depth:
            return
        for event_type, mut_type, base_qual in event_list:
            if event_type == 'good':
                good_kmers[base_qual][mut_type][kmer] += 1
            elif event_type == 'bad':
                bad_kmers[base_qual][mut_type][kmer] += 1

def get_good_and_bad_kmers_w_filter(bam_file, twobit_file, filter_bam_file, basequals, chrom, start, end, min_depth, max_depth, radius):
    finder = MutationFinderWFilter(bam_file, twobit_file, filter_bam_file, basequals)
    good_kmers, bad_kmers = finder.find_mutations(chrom, start, end, min_depth=min_depth, max_depth=max_depth, radius=radius)
    return good_kmers, bad_kmers
