import py2bit
from betterbasequals.utils import reverse_complement, Read, zip_pileups, open_bam_w_index
from betterbasequals.pilup_handlers import get_pileup_count, get_alleles_w_quals, get_alleles_w_corrected_quals


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
        if filter_bam_file is None:
            self.filter_bam_file = None
        else:
            self.filter_bam_file = open_bam_w_index(filter_bam_file)
        self.validation_bam_file = open_bam_w_index(validation_bam_file)

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


    def call_mutations(self, chrom, start, stop):
        if self.filter_bam_file is None:
            self.call_mutations_no_filter(chrom, start, stop)
        else:
            self.call_mutations_with_filter(chrom, start, stop)

    def call_mutations_with_filter(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual_filter=20, prefix=""):
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
                    
    def call_mutations_no_filter(self, chrom, start, stop, mapq=50, mapq_hifi=40, prefix=""):
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
        Hifi_pileup = self.validation_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq_hifi,
            flag_filter=3848,
        )
        
        for pileupcolumn, hifi_pc in zip_pileups(pileup, Hifi_pileup):
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

            corrected_base_quals, n_ref = get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, self.correction_factor)

            n_alt = sum(len(corrected_base_quals[x]) for x in corrected_base_quals)

            hifi_basequals = get_alleles_w_quals(hifi_pc)
            n_hifi_reads = sum(len(hifi_basequals[x]) for x in hifi_basequals)
            
            for A in [x for x in ['A','C','G','T'] if x != ref]:
                if len(corrected_base_quals[A]) == 0:
                    continue
                
                # Variant quality
                corr_var_qual = sum(x[0] for x in corrected_base_quals[A])
                corr_var_qual2 = sum(x[1] for x in corrected_base_quals[A])
                uncorr_var_qual = sum(x[2] for x in corrected_base_quals[A])
                for corrected_Q, uncorrected_Q, base_type in corrected_base_quals[A]:
                    print(chrom, ref_pos, A, int(corrected_Q), uncorrected_Q, base_type, sum(hifi_basequals[A]), n_hifi_reads)
