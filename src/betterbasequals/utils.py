import csv
import sys
#from itertools import product, count
#from operator import add
#from functools import reduce, partial
from collections import defaultdict, Counter
from math import log10
import pysam
import os
import gzip

change_mtypes = ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T')
mtypes = ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T', 'A->A', 'C->C')


complement = str.maketrans("ATCGN", "TAGCN")

def reverse_complement(s):
    return s.translate(complement)[::-1]

def p2phred(p):
    return -10*log10(p)

def phred2p(Q):
    return 10**(-Q/10)

def empirical_bayes(a, b, x, N):
    return (a+x)/(a+b+N)
    #return p2phred(a+x) - p2phred(a+b+N)

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
def mut_type(from_base, to_base, kmer):
    center = len(kmer)//2
    new_kmer = kmer[:center] + from_base + kmer[center+1:]
    if from_base in ['G','T']:
        return f'{ostrand[from_base]}->{ostrand[to_base]}', reverse_complement(new_kmer)
    else:
        return f'{from_base}->{to_base}', new_kmer

def SW_type(base):
    if base in ['G','T']:
        return ostrand[base]
    else:
        return base


def get_average_coverage(bamfile):
    """Quickly calculates the average coverage of a WGS bamfile using the bam index.
    OBS: Only works for bam not cram!

    Args:
        bamfile (pysam.AlignmentFile): bamfile

    Returns:
        int: average coverage over all cromosomes where at least one reads map.
    """
    read_len = 0
    i = 0
    for read in bamfile.fetch():
        read_len = max(read_len, read.query_length)
        i += 1
        if i > 10:
            break
    n_mapped = 0
    genome_len = 0
    for idx_stats in bamfile.get_index_statistics():
        print(idx_stats)
        if idx_stats.mapped > 0:
            n_mapped += idx_stats.mapped
            genome_len += bamfile.get_reference_length(idx_stats.contig)

    return (n_mapped * read_len)/genome_len



# def regions(bed_file):
#     with open(bed_file) as fp:
#         reader = csv.reader(fp, delimiter="\t")
#         for line in reader:
#             yield line[0], int(line[1]), int(line[2])

# def regions_if_sorted(bed_file):
#     with open(bed_file) as fp:
#         reader = csv.reader(fp, delimiter="\t")
#         for line in reader:
#             yield line[0], int(line[1]), int(line[2])

# def generate_mutation_types(k):
#     if k % 2 == 0:
#         raise ValueError("k must be uneven")
#     types = list()
#     mut_pos = (k - 1) // 2
#     for mer in map(
#         partial(reduce, add),
#         product(*(["ATGC"] * mut_pos + ["TC"] + ["ATGC"] * mut_pos)),
#     ):
#         alts = "ATG"
#         if mer[mut_pos] == "T":
#             alts = "AGC"
#         for alt in alts:
#             after = mer[:mut_pos] + alt + mer[mut_pos + 1 :]
#             mut = ">".join([mer[mut_pos],alt])
#             types.append("_".join([mut, mer]))
#     return sorted(types)

# mutation_types_3_mer = dict(zip(generate_mutation_types(3), count()))

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
        self.mapq = pileup_read.alignment.mapping_quality
        #self.NH = pileup_read.alignment.get_tag("NH")
        self.enddist = min(self.pos, len(pileup_read.alignment.query_sequence)-self.pos)
        # Process cigar stats
        cigar_stats = pileup_read.alignment.get_cigar_stats()[0]
        self.has_indel = sum(cigar_stats[1:4]) != 0
        self.has_clip = sum(cigar_stats[4:6]) != 0
        #self.NM = cigar_stats[8]
        self.NM = cigar_stats[10]
        #print(pileup_read.alignment.get_cigar_stats())
    
    def is_good(self, min_enddist=6, max_mismatch = 2):
        return (self.NM <= max_mismatch and
                not self.has_indel and
                not self.has_clip and
                self.enddist >= min_enddist)


class ReadPair:
    def __init__(self, read1, read2):
        self.length = read1.length
        self.allel = read1.allel
        self.start = read2.start if read1.is_reverse else read1.start
        self.min_base_qual = min(read1.base_qual, read2.base_qual)
        self.max_NM = max(read1.NM, read2.NM)
        self.has_indel = 1 if (read1.has_indel or read2.has_indel) else 0
        self.has_clip = 1 if (read1.has_clip or read2.has_clip) else 0




def zip_pileups(*pileups):
    """Iterate through a number of pilup objects and get the columns that are in all.
    Assumes chromosomes are sorted lexicographically in all files
    Yields:
        list : list of pileup columns objects with identical reference pos.
    """
    pileupcolumns = [None] * len(pileups)
    max_chrom = "0"
    while True:
        try:
            for i in range(len(pileups)):
                pileupcolumns[i]  = next(pileups[i])
            while True:
                chrom = max(x.reference_name for x in pileupcolumns)
                assert chrom >= max_chrom, f"chromosomes should be sorted lexicographically, prev:{max_chrom}, current:{chrom}"
                max_chrom = chrom
                max_pos = max(x.reference_pos for x in pileupcolumns)
                for i in range(len(pileups)):
                    while pileupcolumns[i].reference_name < max_chrom:
                        pileupcolumns[i] = next(pileups[i])
                for i in range(len(pileups)):
                    while pileupcolumns[i].reference_name == max_chrom and pileupcolumns[i].reference_pos < max_pos:
                        pileupcolumns[i] = next(pileups[i])
                if all(x.reference_pos == pileupcolumns[0].reference_pos and x.reference_name == pileupcolumns[0].reference_name for x in pileupcolumns):
                    yield pileupcolumns
                    break
        except StopIteration:
            break

def zip_pileups_single_chrom(*pileups):
    """Iterate throug a number of pilup objects and get the columns that are in all.
    Assumes that all pilups are on same chromosome. 
    Yields:
        list : list of pileup columns objects with identical reference pos.
    """
    pileupcolumns = [None] * len(pileups)
    while True:
        try:
            for i in range(len(pileups)):
                pileupcolumns[i]  = next(pileups[i])
            while True:
                max_pos = max(x.reference_pos for x in pileupcolumns)
                for i in range(len(pileups)):
                    while pileupcolumns[i].reference_pos < max_pos:
                        pileupcolumns[i] = next(pileups[i])
                if all(x.reference_pos == pileupcolumns[0].reference_pos for x in pileupcolumns):
                    yield pileupcolumns
                    break
        except StopIteration:
            break


def read_variant_set(variant_file):
    f = open(variant_file)
    S = set()
    for line in f:
        chrom, pos, allele = line.split()
        S.add((chrom, int(pos), allele))
    return S
    #chrom2set = {}
    #for line in f:
    #    chrom, pos, allele = line.split()
    #    if chrom not in chrom2set:
    #        chrom2set[chrom] = set()
    #        chrom2set[chrom].add((int(pos), allele))
    #return chrom2set


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
        #bqual = int(bqual)
        if bqual not in good_kmers:
            good_kmers[bqual] = {}
            for muttype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T', 'C->C', 'A->A'):
                good_kmers[bqual][muttype] = defaultdict(int)
        good_kmers[bqual][mtype][kmer] = int(count)

    bad_kmers = {}
    for line in opts.input_file_bad:
        bqual, mtype, kmer, count = line.split()
        if count == "0":
            continue
        #bqual = int(bqual)
        if bqual not in bad_kmers:
            bad_kmers[bqual] = {}
            for muttype in ('A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T'):
                bad_kmers[bqual][muttype] = defaultdict(int)
        bad_kmers[bqual][mtype][kmer] = int(count)
    return good_kmers, bad_kmers

# def read_kmer_papas(opts):
#     kmer_papas = {}
#     for line in opts.input_file_kmerpapa:
#         bqual, mtype, pattern, correction_factor = line.split()
#         bqual=int(bqual)
#         if bqual not in kmer_papas:
#             kmer_papas[bqual] = {}
#         if mtype not in kmer_papas[bqual]:
#             kmer_papas[bqual][mtype] = {}
#         for context in matches(pattern):
#             kmer_papas[bqual][mtype][context] = float(correction_factor)
#     return kmer_papas

# def read_kmer_papas_for_test(opts):
#     kmer_papas = {}
#     for line in opts.input_file_kmerpapa:
#         bqual, mtype, pattern, correction_factor = line.split()
#         bqual=int(bqual)
#         if bqual not in kmer_papas:
#             kmer_papas[bqual] = {}
#         if mtype not in kmer_papas[bqual]:
#             kmer_papas[bqual][mtype] = {}
#         kmer_papas[bqual][mtype][pattern] = float(correction_factor)
#     return kmer_papas

def read_kmer_papas(opts):
    kmer_papas = {}
    for line in opts.input_file_kmerpapa:
        bqual, mtype, pattern, alpha, beta = line.split()
        alpha = float(alpha)
        beta = float(beta)
        #bqual=int(bqual)
        if bqual not in kmer_papas:
            kmer_papas[bqual] = {}
        if mtype not in kmer_papas[bqual]:
            kmer_papas[bqual][mtype] = {}
        for context in matches(pattern):
            kmer_papas[bqual][mtype][context] = (alpha, beta)
    return kmer_papas

def read_kmer_papas_for_test(opts):
    kmer_papas = {}
    for line in opts.input_file_kmerpapa:
        bqual, mtype, pattern, alpha, beta = line.split()
        alpha = float(alpha)
        beta = float(beta)
        #bqual=int(bqual)
        if bqual not in kmer_papas:
            kmer_papas[bqual] = {}
        if mtype not in kmer_papas[bqual]:
            kmer_papas[bqual][mtype] = {}
        kmer_papas[bqual][mtype][pattern] = p2phred(alpha/(alpha+beta))
    return kmer_papas


def print_kmer_papas(opts, kmer_papas):
    if not opts.output_file_kmerpapa is None:
        for bqual in kmer_papas:
            for mtype in kmer_papas[bqual]:
                for pat in kmer_papas[bqual][mtype]:
                    Q = kmer_papas[bqual][mtype][pat]
                    print(bqual, mtype, pat, Q, file=opts.output_file_kmerpapa)
        opts.output_file_kmerpapa.close()

# def print_good_and_bad(opts, good_kmers, bad_kmers):
#     if not opts.output_file_good is None:
#         for bqual in good_kmers:
#             for mtype in mtypes:
#                 super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
#                 for kmer in matches(super_pattern):
#                     print(bqual, mtype, kmer, good_kmers[bqual][mtype][kmer] , file = opts.output_file_good)
#         opts.output_file_good.close()
#     if not opts.output_file_bad is None:
#         for bqual in bad_kmers:
#             for mtype in mtypes:
#                 super_pattern = 'N'*opts.radius + mtype[0] + 'N'*opts.radius
#                 for kmer in matches(super_pattern):
#                     print(bqual, mtype, kmer, bad_kmers[bqual][mtype][kmer], file = opts.output_file_bad)
#         opts.output_file_bad.close()

def print_good_and_bad(opts, good_kmers, bad_kmers, single_kmers):
    if not opts.output_file_good is None:
        for tup, count in good_kmers.items():
            bqual, mtype, kmer = tup
            print(bqual, mtype, kmer, count , file = opts.output_file_good)
        opts.output_file_good.close()
    if not opts.output_file_bad is None:
        for tup, count in bad_kmers.items():
            bqual, mtype, kmer = tup
            print(bqual, mtype, kmer, count , file = opts.output_file_bad)
        opts.output_file_bad.close()
    if not opts.output_file_single is None:
        for tup, count in single_kmers.items():
            bqual, mtype, kmer = tup
            print(bqual, mtype, kmer, count , file = opts.output_file_single)
        opts.output_file_single.close()

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


class VcfAfReader:
    """A simple VCF-file reader that only parses and returns AF field
    for each position
    """ 

    def __init__(self, vcfname, ignore_indels=True, field_name="AF"):
        # TODO: check if name ends in gz and then open using gzip
        if vcfname.endswith('.gz') or vcfname.endswith('.bgz') or vcfname.endswith('.zip'):
            self.f = gzip.open(vcfname, 'rt')
        else:
            self.f = open(vcfname)
        
        # Read lines until we are past comment lines
        vcf_line = self.f.readline()
        while vcf_line.startswith("#"):
            vcf_line = self.f.readline()
        
        vcfL = vcf_line.split()
        self.vcf_chrom = vcfL[0]
        self.vcf_pos = int(vcfL[1])
        self.vcf_ref = vcfL[3]
        self.vcf_alt = vcfL[4]
        self.vcf_info = vcfL[7]
        self.ignore_indels = ignore_indels
        self.field_name = field_name
        self.last_pos = -1
        self.last_chrom = "0"

    def parse_line(self, vcf_line):
        vcfL = vcf_line.split()
        if len(vcfL) == 0:
            self.vcf_chrom = 'z' # should sort after all posible query chrom names
            return

        if vcfL[0] == self.vcf_chrom:
            assert int(vcfL[1]) >= self.vcf_pos, "ERROR: vcffile positions not sorted!"
        else:
            assert vcfL[0] >= self.vcf_chrom, "ERROR: vcffile chromosomes not sorted!"

        self.vcf_chrom = vcfL[0]
        self.vcf_pos = int(vcfL[1])
        self.vcf_ref = vcfL[3]
        self.vcf_alt = vcfL[4]
        self.vcf_info = vcfL[7]

    def parse_AF(self):
        fn = self.field_name + '='
        L = self.vcf_info.split(";")
        for x in L:
            if x.startswith(fn):
                return float(x[len(fn):])

    def query(self, chrom, pos, ref):
        #Check that current query is after the last query
        assert chrom > self.last_chrom or (chrom == self.last_chrom and pos > self.last_pos), "ERROR: queries should be in sorted order! Please sort bam file."
        self.last_chrom = chrom
        self.last_pos = pos

        while self.vcf_chrom < chrom or (self.vcf_chrom == chrom and self.vcf_pos < pos):
            self.parse_line(self.f.readline())
       
        if (chrom == self.vcf_chrom) and ( self.vcf_pos == pos):
            af_dict = defaultdict(float)
            while((chrom == self.vcf_chrom) and ( self.vcf_pos == pos)):

                #Check that reference matches
                assert self.vcf_ref[0] == ref, f'''
                    ERROR: ref allele in vcffile ({self.vcf_ref[0]}) doesn't match ref({ref}).
                    line: {self.vcf_chrom} {self.vcf_pos} {self.vcf_ref} {self.vcf_alt} {self.vcf_info}.
                    Check genome build!
                    '''

                #ignore indels if ignore_indels option is set:
                if not self.ignore_indels or len(self.vcf_ref)==1 and len(self.vcf_alt)==1:
                    af_dict[self.vcf_alt] = self.parse_AF()

                self.parse_line(self.f.readline())

            return af_dict
        else:
            return None


def print_base_probs(base_probs):
    for A in base_probs:
        print(f"P_{A}->X", f"P_R->X", "count", sep='\t')
        D = Counter()
        for x in base_probs[A]:
            D[x] += 1
        for x in D:
            p_A2X, p_R2X, MQ = x
            print(phred2p(p_A2X),phred2p(p_R2X), D[x], sep = '\t')