from math import log10
import pysam
import py2bit
from betterbasequals.utils import reverse_complement, Read, zip_pileups, open_bam_w_index
from betterbasequals.pilup_handlers import get_pileup_count
from collections import defaultdict


def get_mut_type(ref, papa_ref, alt):
    if ref != papa_ref:
        mtype = papa_ref + '->' + reverse_complement(alt)
    else:
        mtype = papa_ref + '->' + alt
    return mtype



def get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, correction_factor):
    """
    Returns a dictionary that maps from alleles to a list of tuples of the form: (adjusted_BQ, old_BQ, type)
    Where type is:
    1) Match in overlap
    2) Mismatch in overlap
    3) Base with no overlap.
    """

    reads_mem = {}
    base_quals = {'A':[], 'C':[], 'G':[], 'T':[]}

    #tmp: counting ref... can be removed later
    n_ref = 0

    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
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
            ## Have now done so to test --- Can maybe delete later
            if read.allel == mem_read.allel:
                mut_type = get_mut_type(ref, papa_ref, read.allel)

                #Are ignoring ref alleles pt... should adjust later
                if read.allel == ref:
                    n_ref +=1
                    continue
                adjusted_base_qual1 = correction_factor[read.base_qual][mut_type][kmer]
                adjusted_base_qual2 = correction_factor[mem_read.base_qual][mut_type][kmer]
                adjusted_base_qual = adjusted_base_qual1 + adjusted_base_qual2
                unadjusted_base_qual = max(read.base_qual, mem_read.base_qual)
                base_quals[read.allel].append((adjusted_base_qual1 + adjusted_base_qual2, unadjusted_base_qual, 1))
            else:
                #Are ignoring ref alleles pt... should adjust later
                if read.allel != ref:
                    base_quals[read.allel].append((0, read.base_qual, 2))
                else:
                    n_ref += 1
                #Are ignoring ref alleles pt... should adjust later
                if mem_read.allel != ref:
                    base_quals[read.allel].append((0, mem_read.base_qual, 2))
                else:
                    n_ref += 1

        else:            
            reads_mem[read.query_name] = read

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        if read.allel != ref:
            mut_type = get_mut_type(ref, papa_ref, read.allel)
            adjusted_base_qual = correction_factor[read.base_qual][mut_type][kmer]
            base_quals[read.allel].append((adjusted_base_qual, read.base_qual, 3))
        else:
            n_ref += 1
    return base_quals, n_ref


def get_adjustments(pileupcolumn, ref, papa_ref, kmer, correction_factor, change_dict):
    reads_mem = {}

    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
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

            # overlap matches
            if read.allel == mem_read.allel:
                mut_type = get_mut_type(ref, papa_ref, read.allel)

                #Are ignoring ref alleles pt... should adjust later
                if read.allel == ref:
                    continue

                adjusted_base_qual1 = read.base_qual + correction_factor[mut_type][kmer]
                adjusted_base_qual2 = mem_read.base_qual + correction_factor[mut_type][kmer]
                adjusted_base_qual = adjusted_base_qual1 + adjusted_base_qual2
                change_dict[(read.query_name, read.isR1)].append((read.pos, read.base_qual, read.allel, int(adjusted_base_qual1), 1))
                change_dict[(mem_read.query_name, mem_read.isR1)].append((mem_read.pos, mem_read.base_qual, mem_read.allel, int(adjusted_base_qual2), 1))

            else: # overlap mismatch

                #Are ignoring ref alleles pt... could adjust later
                if read.allel != ref:
                    change_dict[(read.query_name, read.isR1)].append((read.pos, read.base_qual, read.allel, 0, 2))
                #Are ignoring ref alleles pt... could adjust later
                if mem_read.allel != ref:
                    change_dict[(mem_read.query_name, mem_read.isR1)].append((mem_read.pos, mem_read.base_qual, mem_read.allel, 0, 2))
        else:            
            reads_mem[read.query_name] = read

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        if read.allel != ref:
            mut_type = get_mut_type(ref, papa_ref, read.allel)
            adjusted_base_qual = read.base_qual + correction_factor[mut_type][kmer]
            change_dict[(read.query_name, read.isR1)].append((read.pos, read.base_qual, read.allel, int(adjusted_base_qual), 3))



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
        bam_file,
        filter_bam_file,
        twobit_file,
        kmer_papa
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)
        self.filter_bam_file = open_bam_w_index(filter_bam_file)

        self.tb = py2bit.open(twobit_file)


        #assumes that the kmer papa has been turned into phred scaled correction factor
        self.correction_factor = kmer_papa

        # Create correction factor dict:

        # self.correction_factor = {mtype:{} for mtype in kmer_papa}
        # print(self.correction)
        # for mtype in kmer_papa:
        #     for kmer in kmer_papa[mtype]:
        #         p = kmer_papa[mtype][kmer]
        #         self.correction_factor[mtype][kmer] = -10*log10(p/(p-1))
        #         #self.correction_factor[mtype][kmer] = -10*log10(p)
                

    def __del__(self):
        self.tb.close()

    def call_mutations(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual_filter=20, radius=3, prefix=""):
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
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name            
            #if not self.bed_query_func(chrom, ref_pos):
            #    continue

            kmer = self.tb.sequence(prefix + chrom, ref_pos- radius, ref_pos + radius + 1)
            ref = kmer[radius]
            if 'N' in kmer or len(kmer)!= 2*radius +1:
                continue            

            if kmer[radius] in ['T', 'G']:
                kmer = reverse_complement(kmer)
                papa_ref = reverse_complement(ref)
                should_reverse = True
            else:
                papa_ref = ref
                should_reverse = False

            assert len(ref) == 1
            if ref not in "ATGC":
                continue

            n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
            N_filter = n_ref_filter + sum(n_alt_filter.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            if N_filter < 25 or N_filter > 55:
                continue

            corrected_base_quals, n_ref = get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, self.correction_factor)

            n_alt = sum(len(corrected_base_quals[x]) for x in corrected_base_quals)
            for A in [x for x in ['A','C','G','T'] if x != ref]:
                if len(corrected_base_quals[A]) == 0:
                    continue
                if n_alt_filter[A] > 0:
                    continue
                # Variant quality
                var_qual = sum(corrected_base_quals[A])
                yield (chrom, ref_pos, ref, A, corrected_base_quals[A], n_alt + n_ref)
                #print(f'{chrom}\t{ref_pos}\t{ref}\t{A}\t{var_qual}\t{len(base_quals)}\t{N}', file=self.outfile)

class MutationValidator:
    def __init__(
        self,
        bam_file,
        filter_bam_file,
        validation_bam_file,
        twobit_file,
        kmer_papa
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)
        self.filter_bam_file = open_bam_w_index(filter_bam_file)

        if not validation_bam_file is None:
            self.validation_bam_file = open_bam_w_index(validation_bam_file)
        else:
            self.validation_bam_file = None

        self.tb = py2bit.open(twobit_file)
        bquals = list(kmer_papa.keys())
        self.radius = len(next(iter(kmer_papa[bquals[0]]["A->C"].keys())))//2

        #assumes that the kmer papa has been turned into phred scaled correction factor
        self.correction_factor = kmer_papa

        # Create correction factor dict:
        # self.correction_factor = {}
        # for mtype in kmer_papa:
        #     for kmer in kmer_papa[mtype]:
        #         if self.radius is None:
        #             self.radius == len(kmer)//2
        #         p = kmer_papa[mtype][kmer]
        #         self.correction_factor[mtype][kmer] = -10*log10(p/(1-p))

    def __del__(self):
        self.tb.close()

    def call_mutations(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual_filter=20, prefix=""):
        pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
            min_base_quality = 1,
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
        Hifi_pileup = self.validation_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq_filter,
            flag_filter=3848,
        )
        
        for pileupcolumn, filter_pc, hifi_pc in zip_pileups(pileup, filter_pileup, Hifi_pileup):
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name            
            #if not self.bed_query_func(chrom, ref_pos):
            #    continue

            kmer = self.tb.sequence(prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
            ref = kmer[self.radius]
            if 'N' in kmer or len(kmer)!= 2*self.radius +1:
                continue            

            if ref in ['T', 'G']:
                kmer = reverse_complement(kmer)
                papa_ref = reverse_complement(ref)
            else:
                papa_ref = ref

            assert len(ref) == 1
            if ref not in "ATGC":
                continue

            n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
            N_filter = n_ref_filter + sum(n_alt_filter.values())
            #TODO: replace hardcoded numbers with values relative to mean coverage
            if N_filter < 25 or N_filter > 55:
                continue

            corrected_base_quals, n_ref = get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, self.correction_factor)

            n_alt = sum(len(corrected_base_quals[x]) for x in corrected_base_quals)

            hifi_basequals = get_alleles_w_quals(hifi_pc)
            n_hifi_reads = sum(len(hifi_basequals[x]) for x in hifi_basequals)
            
            for A in [x for x in ['A','C','G','T'] if x != ref]:
                if len(corrected_base_quals[A]) == 0:
                    continue
                #if n_alt_filter[A] > 0:
                #    continue
                # Variant quality
                corr_var_qual = sum(x[0] for x in corrected_base_quals[A])
                corr_var_qual2 = sum(x[1] for x in corrected_base_quals[A])
                uncorr_var_qual = sum(x[2] for x in corrected_base_quals[A])
                for corrected_Q, uncorrected_Q, base_type in corrected_base_quals[A]:
                    print(chrom, ref_pos, A, int(corrected_Q), uncorrected_Q, base_type, n_alt_filter[A], sum(hifi_basequals[A]), n_hifi_reads)


class BaseAdjuster:
    def __init__(
        self,
        bam_file,
        twobit_file,
        kmer_papa,
        output_bam
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)

        #Open outputfile for writing
        self.out_file = pysam.AlignmentFile(output_bam, "wb", template=self.bam_file)

        self.tb = py2bit.open(twobit_file)

        #assumes that the kmer papa has been turned into phred scaled correction factor
        self.correction_factor = kmer_papa
                

    def __del__(self):
        self.tb.close()

    def call_mutations(self, chrom, start, stop, mapq=50, radius=3, prefix=""):
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

        change_dict = defaultdict(list)

        # Fill dict with changes:
        for pileupcolumn in pileup:
            ref_pos = pileupcolumn.reference_pos
            chrom = pileupcolumn.reference_name            

            kmer = self.tb.sequence(prefix + chrom, ref_pos- radius, ref_pos + radius + 1)
            ref = kmer[radius]
            if 'N' in kmer or len(kmer)!= 2*radius +1:
                continue            

            if kmer[radius] in ['T', 'G']:
                kmer = reverse_complement(kmer)
                papa_ref = reverse_complement(ref)
            else:
                papa_ref = ref

            assert len(ref) == 1
            if ref not in "ATGC":
                continue
            
            get_adjustments(pileupcolumn, ref, papa_ref, kmer, self.correction_factor, change_dict)

        n_corrected_reads = 0
        n_uncorrected_reads = 0
        n_corrections = 0
        n_filtered = 0
        for read in self.bam_file.fetch(chrom, start, stop):
            if read.is_secondary or read.is_supplementary:
                n_filtered += 1
                continue
            read_id = (read.query_name, read.is_read1)
            if read_id in change_dict:
                n_corrected_reads += 1
                for pos, basequal, allele, adjusted_basequal, atype in change_dict[read_id]:
                    n_corrections += 1
                    assert(read.query_sequence[pos] == allele)
                    assert(read.query_qualities[pos] == basequal)
                    read.query_qualities[pos] = adjusted_basequal
            else:
                n_uncorrected_reads += 1

            self.out_file.write(read)
        return n_corrections, n_corrected_reads, n_uncorrected_reads, n_filtered




            