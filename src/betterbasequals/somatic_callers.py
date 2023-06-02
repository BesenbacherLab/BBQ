import py2bit
from scipy.optimize import minimize_scalar
from betterbasequals.utils import eprint, zip_pileups_single_chrom, open_bam_w_index, phred2p, p2phred
from betterbasequals.pilup_handlers import get_alleles_w_probabities_update
from collections import Counter
import scipy.stats

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
        return sum(p2phred((1-p_map_error)*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) +p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)

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
        double_adjustment,
        mapq = 50,
        min_base_qual = 1,
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
        self.double_adjustment = double_adjustment
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix
        self.filter = False

        #TODO: filter_depth variable should be set based on average coverage in filter file.
        self.min_filter_depth = 10
        self.max_filter_depth = 1000000
        self.min_filter_count = 2
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
                filter_alleles = [x for x,y in n_alleles.items() if y >= self.min_filter_count]
                n_calls += self.handle_pileup(pileupcolumn, filter_alleles)
        else:
            for pileupcolumn in pileup:
                n_calls += self.handle_pileup(pileupcolumn, [])
        return n_calls


    def handle_pileup(self, pileupcolumn, filter_alleles):
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
        elif self.method == 'LR':
            base_probs, BQs, n_mismatch, n_double, n_pos, n_neg = \
                get_alleles_w_probabities_update(pileupcolumn, ref, kmer, self.mut_probs, self.prior_N, self.no_update, self.double_adjustment)
            #base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
            for A in base_probs:
                F_list = []
                if A in filter_alleles:
                    continue
                altMQs = [MQ for x,y,MQ in base_probs[A] if x<y]
                if len(altMQs) == 0:
                    continue
                altMQs.sort()
                medianMQ=altMQs[len(altMQs)//2]
                
                if medianMQ < 30:
                    if self.filter:
                        continue
                    else:
                        F_list.append("lowMQ")


                LR, N, N_A, AF = get_LR(base_probs[A])
                # make LR into p-value
                QUAL = int(-LR)
                #if QUAL < 60:
                #    continue
                #p_val = scipy.stats.chi2.sf(-2*LR, 2)
                #if p_val < self.cutoff:
                if QUAL > 0:
                    oldBQ = [x[0] for x in BQs[A]]
                    newBQ = '[' +','.join(f'{x[1]:.1f}' for x in BQs[A]) + ']'
                    
                    oldBQ_str = str(oldBQ)
                    oldBQ.sort()
                    medianBQ=oldBQ[len(oldBQ)//2]
                    
                    n37 = sum(x[0]==37 for x in BQs[A])
                    if medianBQ < 20:
                        if self.filter:
                            continue
                        else:
                            F_list.append("lowBQ")

                    enddist = [x[2] for x in BQs[A]]
                    enddist_str = str(enddist)
                    enddist.sort()
                    median_enddist = enddist[len(enddist)//2]

                    if median_enddist <= 1:
                        if self.filter:
                            continue
                        else:
                            F_list.append("lowEndDist")

                    n_calls += 1
                    if len(F_list) == 0:
                        FILTER = "PASS"
                    else:
                        FILTER = ','.join(F_list)
                    
                    #print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL}\t{FILTER}\tpval={p_val:.3g};LR={LR:.3f};AF={AF:.3g};N={N};N_A={N_A};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={n_mismatch[A]};n_overlap={n_double[A]};MQ={int(medianMQ)}', file=self.outfile)
                    print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL}\t{FILTER}\tAF={AF:.3g};N={N};N_A={N_A};N_A_37={n37};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={n_mismatch[A]};n_overlap={n_double[A]};MQ={int(medianMQ)};alt_strand=[{n_pos[A]},{n_neg[A]}];enddist={enddist_str}', file=self.outfile)

        return n_calls

