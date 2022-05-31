from math import log10
import pysam
import py2bit
from betterbasequals.utils import reverse_complement, Read, zip_pileups
from betterbasequals.pilup_handlers import get_pileup_count
import os

# Mulig PLan:
# zip bam med BQ og bam BAQ så vi har både BQ og BAQ filter for hver allel.

def get_alleles_w_BQ_and_BAQ(pileupcolumn, BAQ_pileupcolumn):

    BAQ_mem = {}
    
    for pileup_read in BAQ_pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

        # fetch read information
        read = Read(pileup_read)
        BAQ_mem[read.query_name] = read.base_qual

    base_quals = {'A':[], 'C':[], 'G':[], 'T':[]}

    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

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
            BQ = read.base_qual
            BAQ = BAQ_mem.pop(read.query_name)    
            base_quals[read.allel].append((BQ, BAQ))

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        adjusted_base_qual = read.base_qual + correction_factor[mtype][kmer]
        base_quals[read.allel].append(adjusted_base_qual)
    return base_quals


def get_alleles_w_quals(pileupcolumn):
    base_quals = {'A':[], 'C':[], 'G':[], 'T':[]}
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

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
        base_quals[read.allel].append(read.base_qual)

    return base_quals


class MutationCaller:
    def __init__(
        self,
        BQ_bam_file,
        BAQ_bam_file,
        filter_bam_file,
        twobit_file,
        outfile,
    ):
        
        self.BQ_bam = open_bam_w_index(BQ_bam_file)
        self.BAQ_bam = open_bam_w_index(BAQ_bam_file)
        self.filter_bam = open_bam_w_index(filter_bam_file)
        self.tb = py2bit.open(twobit_file)
        
        self.outfile = outfile

    def __del__(self):
        self.tb.close()

    def call_mutations(self, chrom, start, stop, min_MQ=40, min_MQ_in_filter=20, min_BQ=90, min_BAQ=40, prefix=""):
        BQ_pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=min_MQ,
            flag_filter=3848,
        )
        BAQ_pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=min_MQ,
            flag_filter=3848,
        )
        filter_pileup = self.filter_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=min_MQ_in_filter,
            flag_filter=3848,
        )
        for BQ_pc, BAQ_pc, filter_pc in zip_pileups(BQ_pileup, BAQ_pileup, filter_pileup):
            ref_pos = BQ_pc.reference_pos
            chrom = BQ_pc.reference_name            
            #if not self.bed_query_func(chrom, ref_pos):
            #    continue

            assert len(ref) == 1
            if ref not in "ATGC":
                continue

            n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
            N_filter = n_ref_filter + sum(n_alt_filter.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            if N_filter < 25 or N_filter > 55:
                continue

            base_quals = get_alleles_w_corrected_quals(pileupcolumn, out_ref, kmer, should_reverse, self.correction_factor)
            N = sum(len(base_quals[x]) for x in base_quals)
            for A in [x for x in ['A','C','G','T'] if x != ref]:
                if len(base_quals[A]) == 0:
                    continue
                if n_alt_filter[A] > 0:
                    continue
                # Variant quality
                var_qual = sum(base_quals)
                print(f'{chrom}\t{ref_pos}\t{ref}\t{A}\t{var_qual}\t{len(base_quals)}\t{N}', file=self.outfile)


def run_mutation_caller(bam_file, filter_bam_file, twobit_file, kmerpapa, outfile, chrom, start, end, radius):
    caller = MutationCaller(bam_file, filter_bam_file, twobit_file, kmer_papa, outfile)
    caller.call_mutations(chrom, start, end, radius=radius)
