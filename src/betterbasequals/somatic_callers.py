import py2bit
from scipy.optimize import minimize_scalar
from betterbasequals.utils import eprint, zip_pileups_single_chrom, open_bam_w_index, phred2p, p2phred, VcfAfReader, mut_type, Read, change_mtypes, SW_type
from collections import Counter, defaultdict
import scipy.stats
import math

def GQ_sum(alt_BQs, ref_BQs):
    return sum(alt_BQs)

def get_LR(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
    # It is ok that I use phred scales insted of the usual natural log.
    # Because the ratio will be the same as log(x)/log(y) == log10(x)/log10(y))
    def p_data_given_mut(alpha):
        return sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x, _ in base_probs)
        
        # If probabilities in base_probs are not phred scaled. It becomes this:
        #return sum(p2phred(alpha * p_a2x + (1-alpha)*p_r2x) for p_a2x, p_r2x in base_probs)
        
        #Maybe we can weigh the read prob using the mapping quality
        #return sum(p2phred((1-p_map_error)*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) +p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)

    res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')
    # Alternative:
    # Integrate out alpha (maybe I should call the getBF (Bayes Factor) instead of getLR then)
    # LL, error = scipy.integrate.quad(p_data_given_mut, 0, 1)
    # Bør nok have prior fordeling på alpha så.
    N = len(base_probs)
    N_A = sum(1 for x,y,_ in base_probs if x<y)
    alpha = res.x

    #if N_A > 1:
    #    print(base_probs)
    #    print([p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs])
    #    eprint(f'alpha={res.x} N={N} N_A={N_A}')
    #    print(res.fun, sum(p_r2x for p_a2x, p_r2x in base_probs),  res.fun - sum(p_r2x for p_a2x, p_r2x in base_probs))
    #LR = res.fun - sum(p2phred((1-p_map_error)*phred2p(p_r2x)+p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)
    LR = res.fun - sum(p_r2x for p_a2x, p_r2x, _ in base_probs)
    return LR, N, N_A, alpha

def get_LR_with_MQ(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
    # It is ok that I use phred scales insted of the usual natural log.
    # Because the ratio will be the same as log(x)/log(y) == log10(x)/log10(y))
    def p_data_given_mut(alpha):
        #return sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs)
        return sum(p2phred((1-phred2p(p_map_error))*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) + phred2p(p_map_error)*0.25) for p_a2x, p_r2x, p_map_error in base_probs)

    res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')

    N = len(base_probs)
    N_A = sum(1 for x,y,_ in base_probs if x<y)
    alpha = res.x

    LR = res.fun - sum(p2phred((1-phred2p(p_map_error))*phred2p(p_r2x)+phred2p(p_map_error)*0.25) for p_a2x, p_r2x, p_map_error in base_probs)
    return LR, N, N_A, alpha


def get_BF(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
    # It is ok that I use phred scales insted of the usual natural log.
    # Because the ratio will be the same as log(x)/log(y) == log10(x)/log10(y))
    def p_data_given_mut(alpha):
        return sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs)
        
        # If probabilities in base_probs are not phred scaled. It becomes this:
        #return sum(p2phred(alpha * p_a2x + (1-alpha)*p_r2x) for p_a2x, p_r2x in base_probs)
        
        #Maybe we can weigh the read prob using the mapping quality
        #return sum(p2phred((1-p_map_error)*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) +p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)

    #res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')
    # Alternative:
    # Integrate out alpha (maybe I should call the getBF (Bayes Factor) instead of getLR then)
    LL, error = scipy.integrate.quad(p_data_given_mut, 0, 1)
    # Bør nok have prior fordeling på alpha så.
    N = len(base_probs)
    N_A = sum(1 for x,y in base_probs if x<y)
    alpha = res.x

    #if N_A > 1:
    #    print(base_probs)
    #    print([p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs])
    #    eprint(f'alpha={res.x} N={N} N_A={N_A}')
    #    print(res.fun, sum(p_r2x for p_a2x, p_r2x in base_probs),  res.fun - sum(p_r2x for p_a2x, p_r2x in base_probs))
    #LR = res.fun - sum(p2phred((1-p_map_error)*phred2p(p_r2x)+p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)
    LR = res.fun - sum(p_r2x for p_a2x, p_r2x in base_probs)
    return LR, N, N_A, alpha



class SomaticMutationCaller:
    def __init__(
        self,
        bam_file,
        filter_bam_file,
        twobit_file,
        kmer_papa,
        outfile,
        method,
        cutoff,
        prior_N,
        no_update,
        N_rate,
        mapq = 50,
        min_base_qual = 1,
        min_filter_count = 2,
        pop_vcf = None,
        min_enddist = 6,
        max_mismatch = 2,
        mean_type = "arithmetric",
        BQ_freq_method = 'global',
        filter_mapq = 20,
        min_base_qual_filter=20, 
        min_depth=1, 
        max_depth=1000000,
        radius=3, 
        prefix="", 
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)
        if filter_bam_file is None:
            self.filter_bam_file = None
        else:
            self.filter_bam_file = open_bam_w_index(filter_bam_file)

        self.tb = py2bit.open(twobit_file)

        # assumes that the kmer papa has been turned into phred scaled probabilities
        self.mut_probs = kmer_papa

        self.outfile = outfile
        self.method = method
        self.BQ_freq_method = BQ_freq_method

        if BQ_freq_method == 'global_by_type':
            self.BQ_freq = {}
            for BQ in kmer_papa:
                self.BQ_freq[BQ] = {}
                for muttype in kmer_papa[BQ]:
                    self.BQ_freq[BQ][muttype] = 0
                    for kmer in kmer_papa[BQ][muttype]:
                        a,b  = kmer_papa[BQ][muttype][kmer]
                        self.BQ_freq[BQ][muttype] += a + b
            
            for muttype in change_mtypes:
                total_sum = sum(self.BQ_freq[BQ][muttype] for BQ in self.BQ_freq)
                for BQ in self.BQ_freq:
                    self.BQ_freq[BQ][muttype] = self.BQ_freq[BQ][muttype]/total_sum

            assert (self.BQ_freq[11]['A->C'] + self.BQ_freq[25]['A->C'] + self.BQ_freq[37]['A->C']) -1.0 < 1e-7
        elif BQ_freq_method == 'global_by_base':
            self.BQ_freq = {}
            for BQ in kmer_papa:
                self.BQ_freq[BQ] = {'A':0, 'C':0}
                for muttype in kmer_papa[BQ]:
                    for kmer in kmer_papa[BQ][muttype]:
                        a,b  = kmer_papa[BQ][muttype][kmer]
                        self.BQ_freq[BQ][muttype[0]] += a + b
            
            for base in ['A', 'C']:
                total_sum = sum(self.BQ_freq[BQ][base] for BQ in self.BQ_freq)
                for BQ in self.BQ_freq:
                    self.BQ_freq[BQ][base] = self.BQ_freq[BQ][base]/total_sum
            assert (self.BQ_freq[11]['A'] + self.BQ_freq[25]['A'] + self.BQ_freq[37]['A']) -1.0 < 1e-7
        elif BQ_freq_method == 'global':
            self.BQ_freq = {}
            for BQ in kmer_papa:
                self.BQ_freq[BQ] = 0
                for muttype in kmer_papa[BQ]:                
                    for kmer in kmer_papa[BQ][muttype]:
                        a,b  = kmer_papa[BQ][muttype][kmer]
                        self.BQ_freq[BQ] += a + b
        
            total_sum = sum(self.BQ_freq.values())
            total_sum += (N_rate*total_sum)/(1-N_rate)
            for BQ in self.BQ_freq:
                self.BQ_freq[BQ] = self.BQ_freq[BQ]/total_sum
            
            assert (self.BQ_freq[11] + self.BQ_freq[25] + self.BQ_freq[37] + N_rate) -1.0 < 1e-7
        else:
            assert False, f'Unknown BQ_freq_method : {BQ_freq_method}'

        #This should be handled in pileup code now.
        #if self.method == 'LR':
            # make sure that all X->X are also in kmerpapa
            # mtype_tups = ('C->C', ('C->A', 'C->G', 'C->T')), ('A->A', ('A->C', 'A->G', 'A->T'))
            # for BQ in self.mut_probs:
            #     for stay_type, change_types in mtype_tups:
            #         self.mut_probs[BQ][stay_type] = {}
            #         for kmer in self.mut_probs[BQ][change_types[0]]:
            #             p = 1.0
            #             for change_type in change_types:
            #                 #TODO: should we use log-sum-exp function for numerical stability?
            #                 #p -= phred2p(self.mut_probs[BQ][change_type][kmer])
            #                 alpha,beta = self.mut_probs[BQ][change_type][kmer]
            #                 p -= alpha/(alpha+beta)
            #             self.mut_probs[BQ][stay_type][kmer] = (p, None)

        self.cutoff = cutoff
        self.prior_N = prior_N
        self.no_update = no_update
        self.N_rate = N_rate

        if not pop_vcf is None:
            self.pop_vcf = VcfAfReader(pop_vcf)
        else:
            self.pop_vcf = None

        self.min_enddist = min_enddist
        self.max_mismatch = max_mismatch
        self.mean_type = mean_type
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix
        self.filter_variants = False
        self.filter_reads = True


        #TODO: filter_depth variable should be set based on average coverage in filter file.
        self.min_filter_depth = 10
        self.max_filter_depth = 1000000
        self.min_filter_count = min_filter_count

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

    def call_all_chroms(self):
        n_calls = 0
        for idx_stats in self.bam_file.get_index_statistics():
            if idx_stats.mapped > 0:
                n_calls += self.call_mutations(idx_stats.contig, None, None)
        return n_calls

    def call_mutations(self, chrom, start, stop):
        pileup = self.bam_file.pileup(
            contig = chrom,
            start = start,
            stop = stop,
            truncate = True,
            max_depth = 1000000,
            min_mapping_quality = self.mapq,
            ignore_overlaps = False,
            flag_require = 0,  # No requirements
            flag_filter = 3840,
            min_base_quality = self.min_base_qual,
        )
        n_calls = 0
        if not self.filter_bam_file is None:
            filter_pileup = self.filter_bam_file.pileup(
                contig = chrom,
                start = start,
                stop = stop,
                truncate = True,
                min_mapping_quality = self.filter_mapq,
                ignore_overlaps = False,
                flag_require = 2,  # proper paired
                flag_filter = 3840,
                min_base_quality = self.min_base_qual_filter,
            )
            for pileupcolumn, filter_pc in zip_pileups_single_chrom(pileup, filter_pileup):
                n_alleles = Counter()
                for pread in filter_pc.pileups:
                    pos = pread.query_position
                    if not pos is None:
                        n_alleles[pread.alignment.query_sequence[pos]] += 1
                N_filter = sum(n_alleles.values())
                if N_filter < self.min_filter_depth or N_filter > self.max_filter_depth:
                    continue
                #filter_alleles = [x for x,y in n_alleles.items() if y >= self.min_filter_count]
                n_calls += self.handle_pileup(pileupcolumn, n_alleles)
        else:
            for pileupcolumn in pileup:
                n_calls += self.handle_pileup(pileupcolumn, None)
        return n_calls


    def handle_pileup(self, pileupcolumn, n_alleles):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%100000 == 0:
            eprint(f"{chrom}:{ref_pos}")            
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue

        n_calls = 0

        if ref_pos-self.radius < 0:
            return 0
        kmer = self.tb.sequence(self.prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
        if 'N' in kmer:
            return 0
        if len(kmer) != 2*self.radius + 1:
            return 0
        ref = kmer[self.radius]
        if ref not in "ATGC":
            return 0

        if self.method == 'poisson':
            raise NotImplementedError
            #make new pileup handler
        elif self.method == 'sum':
            raise NotImplementedError
            #sum = sum(from_R for from_A,from_R in base_probs[A] if from_A < from_R)
        elif self.method in ['LR','LR_with_MQ']:
            base_probs, BQs, n_mismatch, n_double, n_pos, n_neg = \
                self.get_alleles_w_probabities_update(pileupcolumn, ref, kmer)
            #base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

            #no_filter_base_probs, no_filter_BQs, no_filter_n_mismatch, no_filter_n_double, no_filter_n_pos, no_filter_n_neg = \
            #    get_alleles_w_probabities_update(pileupcolumn, ref, kmer, self.mut_probs, self.prior_N, self.no_update, self.double_adjustment, filter_reads=False)

            if not self.pop_vcf is None:
                pop_af = self.pop_vcf.query(chrom, ref_pos+1, ref)
            else:
                pop_af = None

            for A in base_probs:
                F_list = []
                
                if not n_alleles is None and n_alleles[A] > self.min_filter_count:
                    continue
                #if A in filter_alleles:
                #    continue
                altMQs = [MQ for x,y,MQ in base_probs[A] if x<y]
                if len(altMQs) == 0:
                    continue
                altMQs.sort()
                medianMQ=altMQs[len(altMQs)//2]
                
                if medianMQ < 30:
                    if self.filter_variants:
                        continue
                    else:
                        F_list.append("lowMQ")
                if self.method == 'LR':
                    LR, N, N_A, AF = get_LR(base_probs[A])
                elif self.method == 'LR_with_MQ':
                    LR, N, N_A, AF = get_LR_with_MQ(base_probs[A])

                QUAL = -LR
                if self.cutoff is None or QUAL >= self.cutoff:
                    oldBQ_str = '[' + ','.join(x[-1] for x in BQs[A]) + ']'
                    newBQ = '[' +','.join(f'{x[1]:.1f}' for x in BQs[A]) + ']'
                    
                    oldBQ = [x[0] for x in BQs[A]]
                    oldBQ.sort()
                    
                    medianBQ=oldBQ[len(oldBQ)//2]
                    
                    n37 = sum(x[0]==37 for x in BQs[A])

                    if medianBQ < 20:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("lowBQ")

                    enddist = [x[2] for x in BQs[A]]
                    enddist_str = str(enddist)
                    enddist.sort()
                    median_enddist = enddist[len(enddist)//2]

                    if median_enddist <= 1:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("lowEndDist")

                    has_indel = [x[3] for x in BQs[A]]
                    has_clip = [x[4] for x in BQs[A]]
                    NM = [x[5] for x in BQs[A]]
                    #NM_str = str(NM)
                    NM.sort()
                    median_NM = NM[len(NM)//2]
                    frac_indel = sum(has_indel)/len(has_indel)
                    frac_clip = sum(has_clip)/len(has_clip)

                    n37_other = sum(x[0]==37 for alt in 'ACGT' for x in BQs[alt] if alt not in [ref,A])
                    #n37_other_nf = sum(x[0]==37 for alt in 'ACGT' for x in no_filter_BQs[alt] if alt not in [ref,A])

                    n_calls += 1
                    if len(F_list) == 0:
                        FILTER = "PASS"
                    else:
                        FILTER = ','.join(F_list)

                    optional_info = ''
                    if not n_alleles is None:
                        optional_info = f';N_A_f={n_alleles[A]}'
                    if not pop_af is None and pop_af[A] > 1e-7:
                        optional_info += f';pop_AF={pop_af[A]}'

                    #print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL}\t{FILTER}\tpval={p_val:.3g};LR={LR:.3f};AF={AF:.3g};N={N};N_A={N_A};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={n_mismatch[A]};n_overlap={n_double[A]};MQ={int(medianMQ)}', file=self.outfile)
                    #print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL}\t{FILTER}\tAF={AF:.3g};N={N};N_A={N_A};N_A_37={n37};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={no_filter_n_mismatch[A]};n_overlap={no_filter_n_double[A]};MQ={int(medianMQ)};alt_strand=[{n_pos[A]},{n_neg[A]}];enddist={enddist_str};NM={median_NM};frac_indel={frac_indel:.3g};frac_clip={frac_clip:.3g};kmer={kmer};n_other={n37_other};no_filter_n_other={n37_other_nf}{optional_info}', file=self.outfile)
                    print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL:.2f}\t{FILTER}\tAF={AF:.3g};N={N};N_A={N_A};N_A_37={n37};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={n_mismatch[A]};n_overlap={n_double[A]};MQ={int(medianMQ)};alt_strand=[{n_pos[A]},{n_neg[A]}];enddist={enddist_str};NM={median_NM};frac_indel={frac_indel:.3g};frac_clip={frac_clip:.3g};kmer={kmer};n_other={n37_other}{optional_info}', file=self.outfile)

        return n_calls
    
    def get_alleles_w_probabities_update(self, pileupcolumn, ref, ref_kmer):
        """
        Returns a dictionary that maps from allele A to a list of tuples with probability of 
        observing Alt allele A given then read and probability of observing ref allele R given 
        the read. I.e.: base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
        The probabilities are given on phred scale.
        We only considder reads where X_read_i == A or X_read_i == R.
        """

        reads_mem = {}
        seen_alt = set()
        n_mismatch = Counter()
        n_double = Counter()
        n_pos = Counter()
        n_neg = Counter()
        events = {'A':[], 'C':[], 'G':[], 'T':[]}
        double_combinations = set()
        R = ref
        for pileup_read in pileupcolumn.pileups:
            # test for deletion at pileup
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            #TODO: should consider what the right solution is if there is deletion at overlap

            # fetch read information
            read = Read(pileup_read)

            #if filter_reads and not read.is_good():
            #    continue

            # test if read is okay
            if (
                read.allel not in "ATGC"
                or read.start is None
                or read.end is None
            ):
                continue


            # Look for read partner
            if read.query_name in reads_mem:
                # found partner process read pair
                mem_read = reads_mem.pop(read.query_name)

                # We do not trust mismathces in overlaps so we only add to events list in case of match
                if read.allel == mem_read.allel:
                    X = read.allel

                    if self.filter_reads:
                        if (not read.is_good(self.min_enddist, self.max_mismatch)) and mem_read.is_good(self.min_enddist, self.max_mismatch):
                            #considder mem_read single read.
                            reads_mem[read.query_name] = mem_read
                            continue
                        elif (not mem_read.is_good(self.min_enddist, self.max_mismatch)) and read.is_good(self.min_enddist, self.max_mismatch):
                            #considder read single read.
                            reads_mem[read.query_name] = read
                            continue
                        elif (not read.is_good(self.min_enddist, self.max_mismatch)) and (not mem_read.is_good(self.min_enddist, self.max_mismatch)):
                            continue

                    if X == R:
                        alts = [A for A in ['A','C','G','T'] if A!=R]
                    else:
                        alts = [X]
                        seen_alt.add(X)
                        n_pos[X] += 1
                        n_neg[X] += 1

                    for A in alts:
                        #if not no_update:   
                        n_double[A] += 1

                        read_MQ = (read.mapq + mem_read.mapq)/2
                        #if overlap_type == "double":
                        
                        if read.base_qual > mem_read.base_qual:
                            read_BQ = (read.base_qual, mem_read.base_qual)
                            BQ_pair = f'({read.base_qual},{mem_read.base_qual})'
                        else:
                            read_BQ = (mem_read.base_qual, read.base_qual)
                            BQ_pair = f'({mem_read.base_qual},{read.base_qual})'
                        double_combinations.add(read_BQ)
                        enddist = max(read.enddist, mem_read.enddist)
                        has_indel = max(read.has_indel, mem_read.has_indel)
                        has_clip = max(read.has_clip, mem_read.has_clip)
                        NM = max(read.NM, mem_read.NM)
                        #(read.base_qual, mem_read.base_qual)
                        events[A].append((X, read_BQ, read_MQ, enddist, has_indel, has_clip, NM, BQ_pair))

                else: # Mismatch
                    #if not no_update:
                    # TODO: Would it make sense to also count ref_mismatches so that we could
                    # do bayesian update of A->R error rates and not only R->A error rates?
                    if read.allel != ref and mem_read.allel == ref:
                        n_mismatch[read.allel] += 1
                        for A in ['A','C','G','T']:
                            if A == ref:
                                continue
                            n_double[A] += 1
                    if mem_read.allel != ref and read.allel == ref:
                        n_mismatch[mem_read.allel] += 1
                        for A in ['A','C','G','T']:
                            if A == ref:
                                continue
                            n_double[A] += 1
            else:            
                reads_mem[read.query_name] = read

        # Handle reads without partner (ie. no overlap)
        for read in reads_mem.values():
            X = read.allel
            if self.filter_reads and not read.is_good(self.min_enddist, self.max_mismatch):
                continue
                
            if X == R:
                alts = [A for A in ['A','C','G','T'] if A!=R]
            else:

                alts = [X]
                seen_alt.add(X)
                if read.is_reverse:
                    n_neg[X] += 1
                else:
                    n_pos[X] += 1
            
            for A in alts:
                events[A].append((X, read.base_qual, read.mapq, read.enddist, read.has_indel, read.has_clip, read.NM, str(read.base_qual)))
        
        if len(seen_alt) == 0:
            return {}, {}, n_mismatch, n_double, n_pos, n_neg

        new_mut_probs = defaultdict(dict)

        #I only need to calculate probabilities of changing bases from one of the seen alleles.
        # I have to considder change to all bases to calculate stay types (X->X) correctly.
        relevant_bases = [ref] + list(seen_alt)

        # calculate error probabilities for single reads:
        for BQ in self.mut_probs:
            new_mut_probs[BQ] = defaultdict(dict)

            for from_base in relevant_bases:
                p_rest = 1.0
                stay_type, stay_kmer = mut_type(from_base, from_base, ref_kmer)
                for to_base in ['A', 'C', 'G', 'T']:
                    if to_base == from_base:
                        continue
                    change_type, change_kmer = mut_type(from_base, to_base, ref_kmer)
                    alpha, beta = self.mut_probs[BQ][change_type][change_kmer]
                    p_prior = alpha / (alpha + beta)
                    new_p_prior = 0
                    for other_BQ in self.mut_probs:
                        other_alpha, other_beta = self.mut_probs[other_BQ][change_type][change_kmer]
                        other_p_prior = other_alpha / (other_alpha + other_beta)

                        if self.mean_type == "geometric":
                            if self.BQ_freq_method == 'global':
                                new_p_prior += self.BQ_freq[other_BQ] * math.sqrt(other_p_prior*p_prior)
                            elif self.BQ_freq_method == 'global_by_type':
                                new_p_prior += self.BQ_freq[other_BQ][change_type] * math.sqrt(other_p_prior*p_prior)
                            elif self.BQ_freq_method == 'global_by_base':
                                new_p_prior += self.BQ_freq[other_BQ][SW_type(to_base)] * math.sqrt(other_p_prior*p_prior)
                        elif self.mean_type == "arithmetric":
                            if self.BQ_freq_method == 'global':
                                new_p_prior += self.BQ_freq[other_BQ] * ((other_p_prior+p_prior)/2)
                            elif self.BQ_freq_method == 'global_by_type':
                                new_p_prior += self.BQ_freq[other_BQ][change_type] * ((other_p_prior+p_prior)/2)
                            elif self.BQ_freq_method == 'global_by_base':
                                new_p_prior += self.BQ_freq[other_BQ][SW_type(to_base)] * ((other_p_prior+p_prior)/2)
                        else:
                            assert False, f'unknown mean_type parameter: {self.mean_type}'
                    if self.mean_type == "geometric" and self.BQ_freq_method == 'global':
                        new_p_prior += self.N_rate * math.sqrt(p_prior)

                    assert 0.0 < new_p_prior < 1.0, f"prior error probability outside bounds: {new_p_prior}"

                    a = new_p_prior * self.prior_N
                    b = self.prior_N - a
                    if self.no_update or from_base != ref:
                        p_posterior = a/(a + b)
                    else:
                        p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                    p_rest -= p_posterior

                    assert 0.0 < p_posterior < 1.0, f"posterior error probability outside bounds: {p_posterior}"
                    #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                    new_mut_probs[BQ][change_type][change_kmer] = p2phred(p_posterior)

                assert 0.0 < p_rest < 1.0, f"no change posterior error probability outside bounds: {p_rest}"
                new_mut_probs[BQ][stay_type][stay_kmer] = p2phred(p_rest)
    

        # calculate error rates for double reads:
        for BQ1, BQ2 in double_combinations:
            new_mut_probs[(BQ1,BQ2)] = defaultdict(dict)

            for from_base in relevant_bases:
                p_rest = 1.0
                stay_type, stay_kmer = mut_type(from_base, from_base, ref_kmer)
                for to_base in ['A', 'C', 'G', 'T']:
                    if to_base == from_base:
                        continue
                    change_type, change_kmer = mut_type(from_base, to_base, ref_kmer)

                    alpha1, beta1 = self.mut_probs[BQ1][change_type][change_kmer]
                    alpha2, beta2 = self.mut_probs[BQ2][change_type][change_kmer]
                
                    p_prior_1 = alpha1 / (alpha1 + beta1)
                    p_prior_2 = alpha2 / (alpha2 + beta2)

                    if self.mean_type == "geometric":
                        new_p_prior = math.sqrt(p_prior_1 * p_prior_2)
                    elif self.mean_type == "arithmetric":
                        new_p_prior = ((p_prior_1 + p_prior_2)/2)
                    else:
                        assert False, f'unknown mean_type parameter: {self.mean_type}'

                    assert 0.0 < new_p_prior < 1.0, f"prior error probability outside bounds: {new_p_prior}"

                    a = new_p_prior * self.prior_N
                    b = self.prior_N - a
                    if self.no_update or from_base != ref:
                        p_posterior = a/(a + b)
                    else:
                        p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                    p_rest -= p_posterior

                    assert 0.0 < p_posterior < 1.0, f"posterior error probability outside bounds: {p_posterior}"
                    #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                    new_mut_probs[(BQ1,BQ2)][change_type][change_kmer] = p2phred(p_posterior)

                assert 0.0 < p_rest < 1.0, f"no change posterior error probability outside bounds: {p_rest}"
                new_mut_probs[(BQ1,BQ2)][stay_type][stay_kmer] = p2phred(p_rest)    


        posterior_base_probs = {}
        BQs = {'A':[], 'C':[], 'G':[], 'T':[]}
        for A in seen_alt:
            posterior_base_probs[A] = []
            for X, read_BQ, read_MQ, enddist, has_indel, has_clip, NM, BQ_pair in events[A]:
                muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
                muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
                posterior_from_A = new_mut_probs[read_BQ][muttype_from_A][kmer_from_A]
                posterior_from_R = new_mut_probs[read_BQ][muttype_from_R][kmer_from_R]

                posterior_base_probs[A].append((posterior_from_A, posterior_from_R, read_MQ))
                if A==X:
                    if type(read_BQ) == tuple:
                        BQ1, BQ2 = read_BQ
                        read_BQ = max(BQ1,BQ2)
                    BQs[A].append((read_BQ, posterior_from_R, enddist, has_indel, has_clip, NM, BQ_pair))
                
        return posterior_base_probs, BQs, n_mismatch, n_double, n_pos, n_neg


