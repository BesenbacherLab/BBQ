import pysam
import py2bit
import os
import sys
from collections import Counter

from betterbasequals.utils import reverse_complement

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Read:
    def __init__(self, pileup_read):
        pos = pileup_read.query_position
        # set attributes
        self.start = pileup_read.alignment.reference_start
        self.end = pileup_read.alignment.reference_end
        self.allel = pileup_read.alignment.query_sequence[pos]
        self.is_reverse = pileup_read.alignment.is_reverse
        self.base_qual = pileup_read.alignment.query_qualities[pos]
        self.query_name = pileup_read.alignment.query_name
        self.length = abs(pileup_read.alignment.template_length)
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
        
def get_pileup_count(pileupcolumn, ref, min_base_qual, blood=False):
    n_ref = 0
    n_alt = Counter()
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

        if (read.base_qual < min_base_qual):
            continue

        if (not blood) and (read.has_clip == 1 or
            read.has_indel == 1 or
            read.NM >3):
            continue
            
        if read.allel == ref:
            n_ref += 1
        else:
            n_alt[read.allel] += 1
    return n_ref, n_alt



def get_pileup_count_double(pileupcolumn, ref, min_base_qual):
    n_ref = 0
    n_alt = Counter()
    reads_mem = dict()
    has_incompatible = {x:False for x in 'ACGT'}
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
            mem_read = reads_mem[read.query_name]

            # Test if readpair is okay
            if (
                read.is_reverse == mem_read.is_reverse
                or read.length != mem_read.length
                or min(read.base_qual, mem_read.base_qual) < min_base_qual
            ):
                continue

            read_pair = ReadPair(read, mem_read)

            #check read_pair filters
            if (read_pair.has_clip == 1 or
                read_pair.has_indel == 1 or
                read_pair.max_NM >3):
                continue
            
            if read.allel != mem_read.allel:
                has_incompatible[read.allel] = True
                has_incompatible[mem_read.allel] = True
                continue
                
            # Count alleles
            if read.allel == ref:
                n_ref += 1
            else:
                n_alt[read.allel] += 1
                
        else:
            # Unseen read partner store read
            reads_mem[read.query_name] = read
    return n_ref, n_alt, has_incompatible


class MutationFinder:
    def __init__(
        self,
        bam_file,
        filter_bam_file,
        twobit_file,
        #output,
    ):
        # Check if index exists, if not create an index file
        bai_filename1 = f"{bam_file}.bai"
        bai_filename2 = bam_file[:-1] + 'i'
        if not os.path.exists(bai_filename1) and not os.path.exists(bai_filename2):
            print(f"No index file found ({bai_filename1}), generating...")
            pysam.index(bam_file)
        filter_bai_filename1 = f"{filter_bam_file}.bai"
        filter_bai_filename2 = filter_bam_file[:-1] + 'i'
        if not os.path.exists(filter_bai_filename1) and not os.path.exists(filter_bai_filename2):
            print(f"No index file found ({filter_bai_filename1}), generating...")
            pysam.index(filter_bam_file)

        # Open files for reading and writing
        self.bam_file = pysam.AlignmentFile(bam_file, "rb")
        self.filter_bam_file = pysam.AlignmentFile(filter_bam_file, "rb")
        self.tb = py2bit.open(twobit_file)

    def __del__(self):
        self.tb.close()

    def find_mutations(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual=30, min_base_qual_filter=20, radius=3):
        good_kmers = Counter()
        bad_kmers = Counter()
        pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
        )
        filter_pileup = self.filter_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq_filter,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
        )
        for pileupcolumn, filter_pc in zip_pileups(pileup, filter_pileup):
            print(pileupcolumn.reference_pos, filter_pc.reference_pos)
            #print("test", pileupcolumn, filter_pc)
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name            #if not self.bed_query_func(chrom, ref_pos):
            #    continue
            ref = self.tb.sequence(chrom, ref_pos, ref_pos + 1)
            assert len(ref) == 1
            if ref not in "ATGC":
                continue
            
            n_ref, n_alt = get_pileup_count(pileupcolumn, ref, min_base_qual)
            N = n_ref + sum(n_alt.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            if N < 100 or N > 220 or N == n_ref:
                continue

            n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
            N_filter = n_ref_filter + sum(n_alt_filter.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            if N_filter < 25 or N_filter > 55:
                continue

            n_ref_double, n_alt_double, has_incomp = get_pileup_count_double(pileupcolumn, ref, min_base_qual)
            N_double = n_ref_double + sum(n_alt_double.values())
            
            for A in [x for x in ['A','C','G','T'] if x != ref]:
                n_A = n_alt[A]
                if n_A == 0:
                    continue
                n_A_filter = min(n_alt_filter[A], 3)
                if n_A_filter > 0:
                    continue
                n_A_double = n_alt_double[A]
                incompatibility = has_incomp[A]

                if n_A_double > 0 and not incompatibility:
                    mut_type = "good"
                elif n_A_double == 0 and incompatibility:
                    mut_type = "bad"
                else:
                    continue
                
                kmer = self.tb.sequence(chrom, ref_pos- radius, ref_pos + radius + 1)
                if 'N' in kmer:
                    continue
                
                if kmer[radius] in ['T', 'G']:
                    kmer = reverse_complement(kmer)
                    
                if mut_type == "good":
                    good_kmers[kmer] += 1
                elif mut_type == "bad":
                    bad_kmers[kmer] += 1
        return good_kmers, bad_kmers

def get_good_and_bad_kmers(bam_file, filter_bam_file, twobit_file, chrom, start, end, radius):
    finder = MutationFinder(bam_file, filter_bam_file, twobit_file)
    good_kmers, bad_kmers = finder.find_mutations(chrom, start, end, radius=radius)
    return good_kmers, bad_kmers
