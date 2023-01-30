import py2bit
from scipy.optimize import minimize_scalar
from betterbasequals.utils import eprint, zip_pileups, open_bam_w_index, phred2p, p2phred
from betterbasequals.pilup_handlers import get_alleles_w_probabities
import scipy.stats

def GQ_sum(alt_BQs, ref_BQs):
    return sum(alt_BQs)

def get_LR(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
    # It is ok that I use phred scales insted of the usual natural log.
    # Because the ratio will be the same as log(x)/log(y) == log10(x)/log10(y))
    def p_data_given_mut(alpha):
        return sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs)
        
        #Maybe we can weigh the read prob using the mapping quality
        #return sum(p2phred((1-p_map_error)*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) +p_map_error*0.25) for p_a2x, p_r2x, p_map_error in base_probs)

    res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')
    # Alternative:
    # Integrate out alpha (maybe I should call the getBF (Bayes Factor) instead of getLR then)
    # LL, error = scipy.integrate.quad(p_data_given_mut, 0, 1)
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
        mapq = 50,
        min_base_qual = 1,
        filter_mapq = 20,
        min_base_qual_filter=20, 
        min_depth=1, 
        max_depth=5000, 
        radius=3, 
        prefix="", 
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)
        self.filter_bam_file = open_bam_w_index(filter_bam_file)

        self.tb = py2bit.open(twobit_file)

        # assumes that the kmer papa has been turned into phred scaled probabilities
        self.mut_probs = kmer_papa

        self.outfile = outfile
        self.method = method

        if self.method == 'LR':
            # make sure that all X->X are also in kmerpapa
            mtype_tups = ('C->C', ('C->A', 'C->G', 'C->T')), ('A->A', ('A->C', 'A->G', 'A->T'))
            for BQ in self.mut_probs:
                for stay_type, change_types in mtype_tups:
                    self.mut_probs[BQ][stay_type] = {}
                    for kmer in self.mut_probs[BQ][change_types[0]]:
                        p = 1.0
                        for change_type in change_types:
                            #TODO: should we use log-sum-exp function for numerical stability?
                            p -= phred2p(self.mut_probs[BQ][change_type][kmer])
                        self.mut_probs[BQ][stay_type][kmer] = p2phred(p)

        self.cutoff = cutoff
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix

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

    def call_mutations(self, chrom, start, stop):
        pileup = self.bam_file.pileup(
            contig = chrom,
            start = start,
            stop = stop,
            truncate = True,
            min_mapping_quality = self.mapq,
            ignore_overlaps = False,
            flag_require = 2,  # proper paired
            flag_filter = 3848,
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
                flag_filter = 3848,
                min_base_quality = self.min_base_qual_filter,
            )
            for pileupcolumn, filter_pc in zip_pileups(pileup, filter_pileup):
                N_filter = 0
                filter_alleles = set()
                for pread in filter_pc.pileups:
                    pos = pread.query_position
                    if not pos is None:
                        filter_alleles.add(pread.alignment.query_sequence[pos])
                    N_filter += 1

                # TODO: Replace hardcoded values
                if N_filter < 25 or N_filter > 55:
                    continue

                n_calls += self.handle_pileup(pileupcolumn, filter_alleles)
        else:
            for pileupcolumn in pileup:
                n_calls += self.handle_pileup(pileupcolumn, [])
        return n_calls


    def handle_pileup(self, pileupcolumn, filter_alleles):
        ref_pos = pileupcolumn.reference_pos
        chrom = pileupcolumn.reference_name
        if ref_pos%10000 == 0:
            eprint(f"{chrom}:{ref_pos}")            
        #if not self.bed_query_func(chrom, ref_pos):
        #    continue

        n_calls = 0

        if ref_pos-self.radius < 0:
            return 
        kmer = self.tb.sequence(self.prefix + chrom, ref_pos- self.radius, ref_pos + self.radius + 1)
        if 'N' in kmer:
            return
        ref = kmer[self.radius]
        if ref not in "ATGC":
            return

        if self.method == 'poisson':
            base_probs, seen_alts = get_alleles_w_probabities(pileupcolumn, ref, kmer, self.mut_probs)
            #base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
            for A in seen_alts:
                if A in filter_alleles:
                    continue
                LR, N, N_A, AF = get_LR(base_probs[A])
                # make LR into p-value
                p_val = scipy.stats.chi2.sf(-2*LR, 2)
                #print(f'pval=={p_val}')
                if p_val < self.cutoff:
                    print(f'{chrom}\t{ref_pos}\t{ref}\t{A}\tpval={p_val};LR={LR};AF={AF};N={N};N_A={N_A}', file=self.outfile)

            raise NotImplementedError
            #make new pileup handler
        elif self.method == 'sum':
            raise NotImplementedError
            #sum = sum(from_R for from_A,from_R in base_probs[A] if from_A < from_R)
        elif self.method == 'LR':
            base_probs, seen_alts = get_alleles_w_probabities(pileupcolumn, ref, kmer, self.mut_probs)
            #base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
            for A in seen_alts:
                if A in filter_alleles:
                    continue
                LR, N, N_A, AF = get_LR(base_probs[A])
                # make LR into p-value
                p_val = scipy.stats.chi2.sf(-2*LR, 2)
                if p_val < self.cutoff:
                    n_calls += 1
                    print(f'{chrom}\t{ref_pos}\t{ref}\t{A}\tpval={p_val};LR={LR};AF={AF};N={N};N_A={N_A}', file=self.outfile)

        return n_calls

