import csv
import sys
from itertools import product, count
from operator import add
from functools import reduce, partial
import pysam
import os

mtypes = ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T')

complement = str.maketrans("ATCGN", "TAGCN")

def reverse_complement(s):
    return s.translate(complement)[::-1]

def regions(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def regions_if_sorted(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def generate_mutation_types(k):
    if k % 2 == 0:
        raise ValueError("k must be uneven")
    types = list()
    mut_pos = (k - 1) // 2
    for mer in map(
        partial(reduce, add),
        product(*(["ATGC"] * mut_pos + ["TC"] + ["ATGC"] * mut_pos)),
    ):
        alts = "ATG"
        if mer[mut_pos] == "T":
            alts = "AGC"
        for alt in alts:
            after = mer[:mut_pos] + alt + mer[mut_pos + 1 :]
            mut = ">".join([mer[mut_pos],alt])
            types.append("_".join([mut, mer]))
    return sorted(types)

mutation_types_3_mer = dict(zip(generate_mutation_types(3), count()))

code = {'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'R':['A', 'G'],
        'Y':['C', 'T'],
        'S':['G', 'C'],
        'W':['A', 'T'],
        'K':['G', 'T'],
        'M':['A', 'C'],
        'B':['C', 'G', 'T'],
        'D':['A', 'G', 'T'],
        'H':['A', 'C', 'T'],
        'V':['A', 'C', 'G'],
        'N':['A', 'C', 'G', 'T']}

def matches(pattern):
    if len(pattern) == 0:
        yield ''
    else:
        for y in matches(pattern[1:]):
            for x in code[pattern[0]]:
                yield x+y


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Read:
    def __init__(self, pileup_read):
        self.pos = pileup_read.query_position

        # set attributes
        self.start = pileup_read.alignment.reference_start
        self.end = pileup_read.alignment.reference_end
        self.allel = pileup_read.alignment.query_sequence[self.pos]
        self.is_reverse = pileup_read.alignment.is_reverse
        self.base_qual = pileup_read.alignment.query_qualities[self.pos]
        self.query_name = pileup_read.alignment.query_name
        self.length = abs(pileup_read.alignment.template_length)
        self.isR1 = pileup_read.alignment.is_read1
        #self.NH = pileup_read.alignment.get_tag("NH")

        # Process cigar stats
        cigar_stats = pileup_read.alignment.get_cigar_stats()[1]
        self.has_indel = sum(cigar_stats[1:4]) != 0
        self.has_clip = sum(cigar_stats[4:6]) != 0
        self.NM = cigar_stats[10]


class ReadPair:
    def __init__(self, read1, read2):
        self.length = read1.length
        self.allel = read1.allel
        self.start = read2.start if read1.is_reverse else read1.start
        self.min_base_qual = min(read1.base_qual, read2.base_qual)
        self.max_NM = max(read1.NM, read2.NM)
        self.has_indel = 1 if (read1.has_indel or read2.has_indel) else 0
        self.has_clip = 1 if (read1.has_clip or read2.has_clip) else 0


# def zip_pileups(p1, p2, p3):
#     while True:
#         try:
#             pileupcolumn1 = next(p1)
#             pileupcolumn2 = next(p2)
#             pileupcolumn3 = next(p3)
#             while(pileupcolumn1.reference_pos < max(pileupcolumn2.reference_pos, pileupcolumn3.reference_pos)):
#                 pileupcolumn1 = next(p1)
#             while(pileupcolumn2.reference_pos < max(pileupcolumn1.reference_pos, pileupcolumn3.reference_pos)):
#                 pileupcolumn2 = next(p2)
#             while(pileupcolumn3.reference_pos < max(pileupcolumn1.reference_pos, pileupcolumn2.reference_pos)):
#                 pileupcolumn3 = next(p3)
#             if pileupcolumn1.reference_pos == pileupcolumn2.reference_pos == pileupcolumn3.reference_pos:
#                 yield (pileupcolumn1, pileupcolumn2, pileupcolumn3)
#         except StopIteration:
#             break

#TODO: Check that this new general function match the commented one above
def zip_pileups(*pileups):
    while True:
        try:
            pileupcolumns = [None] * len(pileups)
            for i in range(len(pileups)):
                pileupcolumns[i]  = next(pileups[i])
            for i in range(len(pileups)):
                max_pos = max(x.reference_pos for x in pileupcolumns)
                while pileupcolumns[i].reference_pos < max_pos:
                    pileupcolumns[i] = next(pileups[i])
            if all(x.reference_pos == pileupcolumns[0].reference_pos for x in pileupcolumns):
                yield pileupcolumns
        except StopIteration:
            break


def open_bam_w_index(bam_file):
    """Check if bam index exists, if not create an index file.
    Then open and return pysam.AlignmentFile.

    Args:
        bam_file (str): file name

    Returns:
        _type_: _description_
    """
    bai_filename1 = f"{bam_file}.bai"
    bai_filename2 = bam_file[:-1] + 'i'
    if not os.path.exists(bai_filename1) and not os.path.exists(bai_filename2):
        print(f"No index file found ({bai_filename1}), generating...")
        pysam.index(bam_file)
    return pysam.AlignmentFile(bam_file, "rb")

def read_kmers(opts):
    if opts.verbosity > 0:
        eprint("Reading good and bad kmers")
    good_kmers = {}
    for line in opts.input_file_good:
        bqual, mtype, kmer, count = line.split()
        bqual = int(bqual)
        if bqual not in good_kmers:
            good_kmers[bqual] = {}
        if mtype not in good_kmers[bqual]:
            good_kmers[bqual][mtype] = {}
        good_kmers[bqual][mtype][kmer] = int(count)

    bad_kmers = {}
    for line in opts.input_file_bad:
        bqual, mtype, kmer, count = line.split()
        bqual = int(bqual)
        if bqual not in bad_kmers:
            bad_kmers[bqual] = {}
        if mtype not in bad_kmers[bqual]:
            bad_kmers[bqual][mtype] = {}
        bad_kmers[bqual][mtype][kmer] = int(count)
    return good_kmers, bad_kmers

def read_kmer_papas(opts):
    kmer_papas = {}
    for line in opts.input_file_kmerpapa:
        bqual, mtype, pattern, correction_factor = line.split()
        bqual=int(bqual)
        if bqual not in kmer_papas:
            kmer_papas[bqual] = {}
        if mtype not in kmer_papas[bqual]:
            kmer_papas[bqual][mtype] = {}
        for context in matches(pattern):
            kmer_papas[bqual][mtype][context] = float(correction_factor)
    return kmer_papas

def print_kmer_papas(opts, kmer_papas):
    if not opts.output_file_kmerpapa is None:
        for bqual in kmer_papas:
            for mtype in kmer_papas[bqual]:
                for pat in kmer_papas[bqual][mtype]:
                    Q = kmer_papas[bqual][mtype][pat]
                    print(bqual, mtype, pat, Q, file=opts.output_file_kmerpapa)
        opts.output_file_kmerpapa.close()

def print_good_and_bad(opts, good_kmers, bad_kmers):
    if not opts.output_file_good is None:
        for bqual in good_kmers:
            for mtype in mtypes:
                super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
                for kmer in matches(super_pattern):
                    print(bqual, mtype, kmer, good_kmers[bqual][mtype][kmer] , file = opts.output_file_good)
        opts.output_file_good.close()
    if not opts.output_file_bad is None:
        for bqual in bad_kmers:
            for mtype in mtypes:
                super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
                for kmer in matches(super_pattern):
                    print(bqual, mtype, kmer, bad_kmers[bqual][mtype][kmer], file = opts.output_file_bad)
        opts.output_file_bad.close()

def parse_opts_region(opts):
    if opts.region is None:
        opts.chrom = None
        opts.start = None
        opts.end = None
    elif ":" in opts.region:
        opts.chrom, end_points = opts.region.split(':')
        opts.start, opts.end = end_points.split('-')
        opts.start = int(opts.start)
        opts.end = int(opts.end)
    else:
        opts.chrom = opts.region
        opts.start = None
        opts.end = None
