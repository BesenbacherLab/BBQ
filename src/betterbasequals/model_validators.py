import py2bit
from betterbasequals.utils import p2phred, eprint, reverse_complement, Read, zip_pileups_single_chrom, open_bam_w_index
from betterbasequals.pilup_handlers import get_pileup_count, get_alleles_w_quals, get_validation_probabities
from collections import Counter

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
        self.radius = len(next(iter(kmer_papa[bquals[0]]["C->T"].keys())))//2
        self.prefix = ''

        #assumes that the kmer papa has been turned into phred scaled correction factor
        self.mut_probs = kmer_papa
        mtype_tups = ('C->C', ('C->A', 'C->G', 'C->T')), ('A->A', ('A->C', 'A->G', 'A->T'))
        for BQ in self.mut_probs:
            for stay_type, change_types in mtype_tups:
                self.mut_probs[BQ][stay_type] = {}
                for kmer in self.mut_probs[BQ][change_types[0]]:
                    p = 1.0
                    for change_type in change_types:
                        #TODO: should we use log-sum-exp function for numerical stability?
                        #p -= phred2p(self.mut_probs[BQ][change_type][kmer])
                        alpha,beta = self.mut_probs[BQ][change_type][kmer]
                        p -= alpha/(alpha+beta)
                    self.mut_probs[BQ][stay_type][kmer] = (p, None)

        # Create correction factor dict:
        # self.correction_factor = {}
        # for mtype in kmer_papa:
        #     for kmer in kmer_papa[mtype]:
        #         if self.radius is None:
        #             self.radius == len(kmer)//2
        #         p = kmer_papa[mtype][kmer]
        #         self.correction_factor[mtype][kmer] = -10*log10(p/(1-p))
        
        self.min_filter_depth = 10 #25
        self.max_filter_depth = 100 #55


    def __del__(self):
        self.tb.close()

    def call_all_chroms(self):
        for idx_stats in self.bam_file.get_index_statistics():
            if idx_stats.mapped > 0:
                self.call_mutations(idx_stats.contig, None, None)


    def call_mutations(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual_filter=20):
        pileup = self.bam_file.pileup(
            contig = chrom,
            start = start,
            stop = stop,
            truncate = True,
            max_depth = 1000000,
            min_mapping_quality = mapq,
            ignore_overlaps = False,
            flag_require = 0,  # No requirements
            flag_filter = 3840,
            min_base_quality = 2,
        )
        Hifi_pileup = self.validation_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq_filter,
            flag_filter=3848,
        )
        n_calls = 0
        if not self.filter_bam_file is None:
            filter_pileup = self.filter_bam_file.pileup(
                contig = chrom,
                start = start,
                stop = stop,
                truncate = True,
                min_mapping_quality = mapq_filter,
                ignore_overlaps = False,
                flag_require = 2,  # proper paired
                flag_filter = 3840,
                min_base_quality = min_base_qual_filter,
            )
            for pileupcolumn, hifi_pc, filter_pc in zip_pileups_single_chrom(pileup, Hifi_pileup, filter_pileup):
                n_alleles = Counter()
                for pread in filter_pc.pileups:
                    pos = pread.query_position
                    if not pos is None:
                        n_alleles[pread.alignment.query_sequence[pos]] += 1
                    N_filter = sum(n_alleles.values())
                if N_filter < self.min_filter_depth or N_filter > self.max_filter_depth:
                    continue
                filter_alleles = [x for x,y in n_alleles.items() if y >= 5]
                self.handle_pileup(pileupcolumn, hifi_pc, filter_alleles)

        else:
            for pileupcolumn, hifi_pc in zip_pileups_single_chrom(pileup, Hifi_pileup):
                self.handle_pileup(pileupcolumn, Hifi_pileup, ['A', 'C', 'G', 'T'])


    def handle_pileup(self, pileupcolumn, hifi_pc, filter_alleles):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%100000 == 0:
            eprint(f"{chrom}:{ref_pos}")            
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue

        if ref_pos-self.radius < 0:
            return 

        kmer = self.tb.sequence(self.prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
        if 'N' in kmer:
            return 
        ref = kmer[self.radius]
        if ref not in "ATGC":
            return 

        if ref not in filter_alleles:
            return 

        hifi_basequals = get_alleles_w_quals(hifi_pc)
        n_hifi_reads = sum(len(hifi_basequals[x]) for x in hifi_basequals)
            
        if n_hifi_reads > 100:
            return

        base_probs, seen_alts, n_mismatch, n_double, n_alt = get_validation_probabities(pileupcolumn, ref, kmer, self.mut_probs)
        
        all_mismatch = sum(n_mismatch.values())
        for A in seen_alts:
            seen_hifi = sum(1 for x in hifi_basequals[A] if x>30) > 4
            for alpha, beta, BQ, muttype, atype in base_probs[A]:
                if n_alt[A]>5:
                    continue
                print(BQ, alpha, beta, muttype, atype, int(seen_hifi), n_mismatch[A], all_mismatch, n_double[A], n_alt[A])

    #             # Variant quality
    #             #corr_var_qual = sum(x[0] for x in corrected_base_quals[A])
    #             #corr_var_qual2 = sum(x[1] for x in corrected_base_quals[A])
    #             #uncorr_var_qual = sum(x[2] for x in corrected_base_quals[A])
    #             n_alt = sum(1 for x in corrected_base_quals[A] if x[1]>35)
    #             seen_hifi = sum(1 for x in hifi_basequals[A] if x>80) > 0
    #             for corrected_Q, uncorrected_Q, base_type in corrected_base_quals[A]:
    #                 print(int(corrected_Q), uncorrected_Q, base_type, seen_hifi, n_mismatch[A], all_mismatch, n_double, n_alt)
    

    # def call_mutations(self, chrom, start, stop):
    #     if self.filter_bam_file is None:
    #         self.call_mutations_no_filter(chrom, start, stop)
    #     else:
    #         self.call_mutations_with_filter(chrom, start, stop)

    # def call_mutations_with_filter(self, chrom, start, stop, mapq=50, mapq_filter=20, min_base_qual_filter=20, prefix=""):
    #     pileup = self.bam_file.pileup(
    #         contig=chrom,
    #         start=start,
    #         stop=stop,
    #         truncate=True,
    #         max_depth = 1000000,
    #         min_mapping_quality=mapq,
    #         ignore_overlaps=False,
    #         flag_require=0,  # proper paired
    #         flag_filter=3848,
    #         min_base_quality = 1,
    #     )
    #     filter_pileup = self.filter_bam_file.pileup(
    #         contig=chrom,
    #         start=start,
    #         stop=stop,
    #         truncate=True,
    #         min_mapping_quality=mapq_filter,
    #         ignore_overlaps=False,
    #         flag_require=2,  # proper paired
    #         flag_filter=3848,
    #     )
    #     Hifi_pileup = self.validation_bam_file.pileup(
    #         contig=chrom,
    #         start=start,
    #         stop=stop,
    #         truncate=True,
    #         min_mapping_quality=mapq_filter,
    #         flag_filter=3848,
    #     )
        
    #     for pileupcolumn, filter_pc, hifi_pc in zip_pileups_single_chrom(pileup, filter_pileup, Hifi_pileup):
    #         ref_pos = pileupcolumn.reference_pos
    #         chrom = pileupcolumn.reference_name            
    #         #if not self.bed_query_func(chrom, ref_pos):
    #         #    continue

    #         kmer = self.tb.sequence(prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
    #         ref = kmer[self.radius]
    #         if 'N' in kmer or len(kmer)!= 2*self.radius +1:
    #             continue            

    #         if ref in ['T', 'G']:
    #             kmer = reverse_complement(kmer)
    #             papa_ref = reverse_complement(ref)
    #         else:
    #             papa_ref = ref

    #         assert len(ref) == 1
    #         if ref not in "ATGC":
    #             continue

    #         n_ref_filter, n_alt_filter = get_pileup_count(filter_pc, ref, min_base_qual_filter, blood=True)
    #         N_filter = n_ref_filter + sum(n_alt_filter.values())
    #         #TODO: replace hardcoded numbers with values relative to mean coverage
    #         if N_filter < 25 or N_filter > 55 or n_ref_filter < 5:
    #             continue

    #         corrected_base_quals, n_ref, n_mismatch, n_double = get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, self.correction_factor)

    #         #n_alt = sum(len(corrected_base_quals[x]) for x in corrected_base_quals)

    #         hifi_basequals = get_alleles_w_quals(hifi_pc)
    #         n_hifi_reads = sum(len(hifi_basequals[x]) for x in hifi_basequals)
            
    #         if n_hifi_reads > 100:
    #             continue

    #         all_mismatch = sum(n_mismatch.values())            

    #         for A in [x for x in ['A','C','G','T'] if x != ref]:    
    #             # Variant quality
    #             #corr_var_qual = sum(x[0] for x in corrected_base_quals[A])
    #             #corr_var_qual2 = sum(x[1] for x in corrected_base_quals[A])
    #             #uncorr_var_qual = sum(x[2] for x in corrected_base_quals[A])
    #             n_alt = sum(1 for x in corrected_base_quals[A] if x[1]>35)
    #             seen_hifi = sum(1 for x in hifi_basequals[A] if x>80) > 0
    #             for corrected_Q, uncorrected_Q, base_type in corrected_base_quals[A]:
    #                 print(int(corrected_Q), uncorrected_Q, base_type, seen_hifi, n_mismatch[A], all_mismatch, n_double, n_alt)
                    
    # def call_mutations_no_filter(self, chrom, start, stop, mapq=50, mapq_hifi=40, prefix=""):
    #     pileup = self.bam_file.pileup(
    #         contig=chrom,
    #         start=start,
    #         stop=stop,
    #         truncate=True,
    #         max_depth = 1000000,
    #         min_mapping_quality=mapq,
    #         ignore_overlaps=False,
    #         flag_require=2,  # proper paired
    #         flag_filter=3848,
    #         min_base_quality = 1,
    #     )
    #     Hifi_pileup = self.validation_bam_file.pileup(
    #         contig=chrom,
    #         start=start,
    #         stop=stop,
    #         truncate=True,
    #         min_mapping_quality=mapq_hifi,
    #         flag_filter=3848,
    #     )
        
    #     for pileupcolumn, hifi_pc in zip_pileups_single_chrom(pileup, Hifi_pileup):
    #         ref_pos = pileupcolumn.reference_pos
    #         chrom = pileupcolumn.reference_name            
    #         #if not self.bed_query_func(chrom, ref_pos):
    #         #    continue

    #         kmer = self.tb.sequence(prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
    #         ref = kmer[self.radius]
    #         if 'N' in kmer or len(kmer)!= 2*self.radius +1:
    #             continue            

    #         if ref in ['T', 'G']:
    #             kmer = reverse_complement(kmer)
    #             papa_ref = reverse_complement(ref)
    #         else:
    #             papa_ref = ref

    #         assert len(ref) == 1
    #         if ref not in "ATGC":
    #             continue

    #         corrected_base_quals, n_ref, n_mismatch, n_double = get_alleles_w_corrected_quals(pileupcolumn, ref, papa_ref, kmer, self.correction_factor)

    #         #n_alt = sum(len(corrected_base_quals[x]) for x in corrected_base_quals)

    #         hifi_basequals = get_alleles_w_quals(hifi_pc)
    #         #n_hifi_reads = sum(len(hifi_basequals[x]) for x in hifi_basequals)
    #         all_mismatch = sum(n_mismatch.values())            

    #         for A in [x for x in ['A','C','G','T'] if x != ref]:    
    #             # Variant quality
    #             #corr_var_qual = sum(x[0] for x in corrected_base_quals[A])
    #             #corr_var_qual2 = sum(x[1] for x in corrected_base_quals[A])
    #             #uncorr_var_qual = sum(x[2] for x in corrected_base_quals[A])

    #             for corrected_Q, uncorrected_Q, base_type in corrected_base_quals[A]:
    #                 print(int(corrected_Q), uncorrected_Q, base_type, sum(1 for x in hifi_basequals[A] if x>80), n_mismatch[A], all_mismatch, n_double)
    #                 #print(chrom, ref_pos, A, int(corrected_Q), uncorrected_Q, base_type, sum(hifi_basequals[A]), n_mismatch, n_hifi_reads)
