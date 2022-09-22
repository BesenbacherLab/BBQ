import pysam
import py2bit
import os
from collections import Counter
from betterbasequals.utils import reverse_complement, mtypes, zip_pileups, eprint, open_bam_w_index
from betterbasequals.pilup_handlers import get_pileup_count, get_pileup_count_double  

class MutationFinderWFilter:
    def __init__(
        self,
        bam_file,
        filter_bam_file,
        twobit_file,
        #output,
    ):
        self.bam_file = open_bam_w_index(bam_file)
        self.filter_bam_file = open_bam_w_index(filter_bam_file)
        self.tb = py2bit.open(twobit_file)

    def __del__(self):
        self.tb.close()

    def find_mutations(self, chrom, start, stop, mapq=50, min_base_qual=30, filter_mapq =20, min_base_qual_filter=20, min_depth=1, max_depth=5000, radius=3, prefix=""):
        good_kmers = {x:Counter() for x in mtypes}
        bad_kmers = {x:Counter() for x in mtypes}
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
            min_mapping_quality=filter_mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
        )
        for pileupcolumn, filter_pc in zip_pileups(pileup, filter_pileup):
            #print(pileupcolumn.reference_pos, filter_pc.reference_pos)
            #print("test", pileupcolumn, filter_pc)
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name
            if ref_pos%10000 ==0:
                eprint(f"{chrom}:{ref_pos}")            
            #if not self.bed_query_func(chrom, ref_pos):
            #    continue
            ref = self.tb.sequence(prefix + chrom, ref_pos, ref_pos + 1)
            assert len(ref) == 1
            if ref not in "ATGC":
                continue

            n_ref, n_alt = get_pileup_count(pileupcolumn, ref, min_base_qual)
            N = n_ref + sum(n_alt.values())
            #TODO: can I learn good filter values from data.
            # fx just have values relative to mean. maybe with poison stderr.
            if N < min_depth or N > max_depth or N == n_ref:
                continue

            n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
            N_filter = n_ref_filter + sum(n_alt_filter.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            #print(N_filter)
            if N_filter < 25 or N_filter > 55:
                continue

            n_ref_double, n_alt_double, has_incomp = get_pileup_count_double(pileupcolumn, ref, min_base_qual)
            N_double = n_ref_double + sum(n_alt_double.values())
            
            #print(N_double)
            #if N_double == 0:
            #    continue
                #print(chrom, ref_pos)

            kmer = self.tb.sequence(prefix + chrom, ref_pos- radius, ref_pos + radius + 1)
            if 'N' in kmer:
                continue
            
            should_reverse = False
            if kmer[radius] in ['T', 'G']:
                kmer = reverse_complement(kmer)
                should_reverse = True

            for A in [x for x in ['A','C','G','T'] if x != ref]:
                n_A = n_alt[A]
                if n_A == 0:
                    continue
                #n_A_filter = min(n_alt_filter[A], 3)
                #if n_A_filter > 0:
                #    continue
                n_A_double = n_alt_double[A]
                incompatibility = has_incomp[A]

                if n_A_double > 0 and not incompatibility:
                    mut_type = "good"
                elif n_A_double == 0 and incompatibility:
                    mut_type = "bad"
                else:
                    continue
                
                if should_reverse:
                    A = reverse_complement(A)

                mtype = kmer[radius] + '->' + A
                assert mtype in mtypes
                if mut_type == "good":
                    good_kmers[mtype][kmer] += 1
                elif mut_type == "bad":
                    bad_kmers[mtype][kmer] += 1
        return good_kmers, bad_kmers

def get_good_and_bad_kmers_w_filter(bam_file, filter_bam_file, twobit_file, chrom, start, end, min_depth, max_depth, radius):
    finder = MutationFinderWFilter(bam_file, filter_bam_file, twobit_file)
    good_kmers, bad_kmers = finder.find_mutations(chrom, start, end, min_depth=min_depth, max_depth=max_depth, radius=radius)
    return good_kmers, bad_kmers


class MutationFinder:
    def __init__(
        self,
        bam_file,
        twobit_file,
    ):
        self.bam_file = open_bam_w_index(bam_file)
        self.tb = py2bit.open(twobit_file)

    def __del__(self):
        self.tb.close()

    def find_mutations(self, chrom, start, stop, mapq=50, min_base_qual=30, min_depth=1, max_depth=5000, radius=3, prefix=""):
        good_kmers = {x:Counter() for x in mtypes}
        bad_kmers = {x:Counter() for x in mtypes}
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
        for pileupcolumn in pileup:
            #print(pileupcolumn.reference_pos, filter_pc.reference_pos)
            #print("test", pileupcolumn, filter_pc)
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name
            if ref_pos%10000 ==0:
                eprint(f"{chrom}:{ref_pos}")            
            #if not self.bed_query_func(chrom, ref_pos):
            #    continue
            ref = self.tb.sequence(prefix + chrom, ref_pos, ref_pos + 1)
            assert len(ref) == 1
            if ref not in "ATGC":
                continue
            
            n_ref, n_alt = get_pileup_count(pileupcolumn, ref, min_base_qual)
            N = n_ref + sum(n_alt.values())
            #TODO: can I learn good filter values from data.
            # fx just have values relative to mean. maybe with poison stderr.
            if N < min_depth or N > max_depth or N == n_ref:
                continue

            n_ref_double, n_alt_double, has_incomp = get_pileup_count_double(pileupcolumn, ref, min_base_qual)
            #N_double = n_ref_double + sum(n_alt_double.values())
            
            kmer = self.tb.sequence(prefix + chrom, ref_pos- radius, ref_pos + radius + 1)
            if 'N' in kmer:
                continue
            
            should_reverse = False
            if kmer[radius] in ['T', 'G']:
                kmer = reverse_complement(kmer)
                should_reverse = True

            for A in [x for x in ['A','C','G','T'] if x != ref]:
                n_A = n_alt[A]
                if n_A == 0:
                    continue
                n_A_double = n_alt_double[A]
                incompatibility = has_incomp[A]

                if n_A_double > 0 and not incompatibility:
                    mut_type = "good"
                elif n_A_double == 0 and incompatibility:
                    mut_type = "bad"
                else:
                    continue
                
                if should_reverse:
                    A = reverse_complement(A)

                mtype = kmer[radius] + '->' + A
                assert mtype in mtypes
                if mut_type == "good":
                    good_kmers[mtype][kmer] += 1
                elif mut_type == "bad":
                    bad_kmers[mtype][kmer] += 1
        return good_kmers, bad_kmers

def get_good_and_bad_kmers(bam_file,  twobit_file, chrom, start, end, min_depth, max_depth, radius):
    finder = MutationFinder(bam_file, twobit_file)
    good_kmers, bad_kmers = finder.find_mutations(chrom, start, end, min_depth=min_depth, max_depth=max_depth, radius=radius)
    return good_kmers, bad_kmers

