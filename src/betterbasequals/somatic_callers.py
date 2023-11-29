import py2bit
from scipy.optimize import minimize_scalar
from betterbasequals.utils import eprint, zip_pileups_single_chrom, open_bam_w_index, phred2p, p2phred, VcfAfReader, mut_type, Read, change_mtypes, SW_type, SortedBed
from collections import Counter, defaultdict
import scipy.stats
import math

def GQ_sum(alt_BQs, ref_BQs):
    return sum(alt_BQs)

def get_LR(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    def p_data_given_mut(alpha):
        return - sum(math.log(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x, _ in base_probs)
        
        # If probabilities in base_probs are not phred scaled. It becomes this:
        #return sum(p2phred(alpha * p_a2x + (1-alpha)*p_r2x) for p_a2x, p_r2x in base_probs)
        
    res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')

    alpha = res.x
    LR = - res.fun + p_data_given_mut(0)
    return LR, alpha


def get_maxLR_with_MQ(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    LL_with_alt = sum(math.log((1-phred2p(p_map_error))*max(phred2p(p_a2x),phred2p(p_r2x)) + phred2p(p_map_error)*0.5) for p_a2x, p_r2x, p_map_error in base_probs if p_r2x > p_a2x)
    LL_no_alt = sum(math.log((1-phred2p(p_map_error))*phred2p(p_r2x) + phred2p(p_map_error)*0.5) for p_a2x, p_r2x, p_map_error in base_probs if p_r2x > p_a2x)

    LR = LL_with_alt - LL_no_alt
    return LR


def get_LR_with_MQ(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    def p_data_given_mut(alpha):
        #return sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x in base_probs)
        return - sum(math.log((1-phred2p(p_map_error))*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) + phred2p(p_map_error)*0.5) for p_a2x, p_r2x, p_map_error in base_probs)

    res = minimize_scalar(p_data_given_mut, bounds=(0, 1), method='bounded')
    alpha = res.x
    LR = - res.fun + p_data_given_mut(0)
    return LR, alpha


def get_BF(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    def p_data_given_mut(alpha):
        return phred2p(sum(p2phred(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) for p_a2x, p_r2x, _ in base_probs))
        
    # Integrate out alpha 
    p_no_alpha, error = scipy.integrate.quad(p_data_given_mut, 0, 1)

    try:
        log_BF = math.log(p_no_alpha) - math.log(p_data_given_mut(0))
    except:
        log_BF = -1000

    return log_BF

def get_BF_with_MQ(base_probs):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    def p_data_given_mut(alpha):
        return phred2p(sum(p2phred((1-phred2p(p_map_error))*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) + phred2p(p_map_error)*0.5) for p_a2x, p_r2x, p_map_error in base_probs))

    # Integrate out alpha
    p_no_alpha, error = scipy.integrate.quad(p_data_given_mut, 0, 1)

    try:
        log_BF = math.log(p_no_alpha) - math.log(p_data_given_mut(0))
    except:
        log_BF = -1000

    return log_BF

def get_BF_with_MQ_and_Prior(base_probs, a=1, b=2):
    #base_probs = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]

    def p_data_given_mut(alpha):
        return scipy.stats.beta.pdf(alpha, a, b) * phred2p(sum(p2phred((1-phred2p(p_map_error))*(alpha * phred2p(p_a2x) + (1-alpha)*phred2p(p_r2x)) + phred2p(p_map_error)*0.5) for p_a2x, p_r2x, p_map_error in base_probs))

    # Integrate out alpha
    p_no_alpha, error = scipy.integrate.quad(p_data_given_mut, 0, 1)

    try:
        log_BF =  math.log(p_no_alpha) - math.log(p_data_given_mut(0))
    except:
        log_BF = - 1000

    return log_BF


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
        mapq = 50,
        min_base_qual = 1,
        min_filter_count = 2,
        pop_vcf = None,
        min_enddist = 6,
        max_mismatch = 2,
        max_NM_diff = 1,
        filter_mapq = 20,
        min_base_qual_filter = 20, 
        min_depth = 1, 
        max_depth = 1000000,
        radius = 3, 
        prefix = "",
        bed_file = None,
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

        self.cutoff = cutoff
        self.prior_N = prior_N
        self.no_update = no_update

        if not pop_vcf is None:
            self.pop_vcf = VcfAfReader(pop_vcf)
        else:
            self.pop_vcf = None

        self.min_enddist = min_enddist
        self.max_mismatch = max_mismatch
        self.max_NM_diff = max_NM_diff
        self.mapq = mapq
        self.min_base_qual = min_base_qual
        self.filter_mapq = filter_mapq
        self.min_base_qual_filter = min_base_qual_filter
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.radius = radius
        self.prefix = prefix
        self.filter_variants = False

        #TODO: filter_depth variable should be set based on average coverage in filter file.
        self.min_filter_depth = 10
        self.max_filter_depth = 1000000
        self.min_filter_count = min_filter_count

        if not bed_file is None:
            self.bed = SortedBed(bed_file)
        else:
            self.bed = None
                

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
            min_mapping_quality = 0, #self.mapq,
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

        if not self.bed is None and self.bed.query(chrom, ref_pos):
            return 0

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
        elif self.method in ['LR','LR_with_MQ', 'BF', 'BF_with_MQ', 'BF_with_MQ_and_Prior', 'maxLR_with_MQ']:
            base_probs, BQs, n_mismatch, n_double, n_pos, n_neg, n_filtered, n_nonfiltered = \
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
                medianMQ = altMQs[len(altMQs)//2]
                
                N_total = sum(n_filtered.values()) + sum(n_nonfiltered.values())

                if medianMQ < 30:
                    if self.filter_variants:
                        continue
                    else:
                        F_list.append("lowMQ")

                refMQs = [MQ for x,y,MQ in base_probs[A] if x>y]
                refMQs.sort()
                
                if len(refMQs)>0:
                    ref_medianMQ = refMQs[len(refMQs)//2]
                else:
                    ref_medianMQ = 0
                    continue

                AF = None
                if self.method == 'LR':
                    QUAL, AF = get_LR(base_probs[A])
                elif self.method == 'LR_with_MQ':
                    QUAL, AF = get_LR_with_MQ(base_probs[A])
                elif self.method == 'maxLR_with_MQ':
                    QUAL = get_maxLR_with_MQ(base_probs[A])
                elif self.method == 'BF':
                    QUAL = get_BF(base_probs[A])
                elif self.method == 'BF_with_MQ':
                    QUAL = get_BF_with_MQ(base_probs[A])
                elif self.method == 'BF_with_MQ_and_Prior':
                    QUAL = get_BF_with_MQ_and_Prior(base_probs[A])

                if self.cutoff is None or QUAL >= self.cutoff:

                    N = len(base_probs[A])
                    N_A = sum(1 for x,y,_ in base_probs[A] if x<y)
                    if AF is None:
                        AF = N_A/N

                    oldBQ_str = '[' + ','.join(x[-1] for x in BQs[A]) + ']'
                    newBQ = '[' +','.join(f'{x[1]:.1f}' for x in BQs[A]) + ']'
                    
                    oldBQ = [x[0] for x in BQs[A]]
                    oldBQ.sort()
                    
                    medianBQ = oldBQ[len(oldBQ)//2]
                    
                    n37 = sum(x[0]==37 for x in BQs[A])
                    a = n_filtered[A]
                    b = n_nonfiltered[A]
                    c = n_filtered[ref]
                    d = n_nonfiltered[ref]
                    if a*d <= b*c:
                        filter_allele_pval = 1.0
                    else:
                        filter_allele_pval = scipy.stats.fisher_exact([[a,b],[c,d]], alternative='greater')[1]
                    filter_allele_table = f'[{a},{b},{c},{d}]'

                    if filter_allele_pval < 0.05:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("FilterBias")

                    indel_reads = n_filtered['-']

                    if medianBQ < 20:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("lowBQ")

                    if min(n_pos[A], n_neg[A]) == 0 and max(n_pos[A], n_neg[A])>4:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("StrandBias")


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
                    NM_alt = [x[5] for x in BQs[A]]
                    NM_ref = [x[5] for x in BQs[ref]]
                    #NM_str = str(NM)
                    NM_alt.sort()
                    NM_ref.sort()
                    median_NM_alt = NM_alt[len(NM_alt)//2]
                    min_NM_alt = min(NM_alt)

                    if len(NM_ref)>0:
                        median_NM_ref = NM_ref[len(NM_ref)//2]
                        min_NM_ref = min(NM_ref)
                        quartile_NM_ref = NM_ref[(len(NM_ref)*3)//4]
                    else:
                        median_NM_ref = 0
                        min_NM_ref = 0
                    
                    #if min_NM_alt - median_NM_ref > self.max_NM_diff:
                    if median_NM_alt - median_NM_ref > self.max_NM_diff:
                        if self.filter_variants:
                            continue
                        else:
                            F_list.append("manyAltMismatches")

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
                    print(f'{chrom}\t{ref_pos+1}\t.\t{ref}\t{A}\t{QUAL:.2f}\t{FILTER}\tAF={AF:.3g};N={N};N_A={N_A};N_A_37={n37};N_total={N_total};oldBQ={oldBQ_str};newBQ={newBQ};n_mismatch={n_mismatch[A]};n_overlap={n_double[A]};MQ={int(medianMQ)};alt_strand=[{n_pos[A]},{n_neg[A]}];enddist={enddist_str};median_alt_NM={median_NM_alt};median_ref_NM={median_NM_ref};min_alt_NM={min_NM_alt};min_ref_NM={min_NM_ref};quartile_NM_ref={quartile_NM_ref};frac_indel={frac_indel:.3g};frac_clip={frac_clip:.3g};kmer={kmer};n_other={n37_other};ref_medianMQ={ref_medianMQ};filter_vs_allele={filter_allele_table};filter_allele_pval={filter_allele_pval:.3g};indel_reads={indel_reads}{optional_info}', file=self.outfile)

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
        events = {'A':{}, 'C':{}, 'G':{}, 'T':{}}

        n_filtered = Counter()
        n_nonfiltered = Counter()

        observed_BQs = set()
        R = ref
        for pileup_read in pileupcolumn.pileups:
           
            # fetch read information
            read = Read(pileup_read)

            # Look for read partner
            if read.query_name in reads_mem:
                # found partner process read pair
                mem_read = reads_mem.pop(read.query_name)

                # We do not trust mismathces in overlaps so we only add to events list in case of match
                if read.allel == mem_read.allel:
                    X = read.allel

                    if (not read.is_good(self.min_enddist, self.max_mismatch, self.mapq)) and mem_read.is_good(self.min_enddist, self.max_mismatch, self.mapq):
                        #considder mem_read single read.
                        reads_mem[read.query_name] = mem_read
                        continue
                    elif (not mem_read.is_good(self.min_enddist, self.max_mismatch, self.mapq)) and read.is_good(self.min_enddist, self.max_mismatch, self.mapq):
                        #considder read single read.
                        reads_mem[read.query_name] = read
                        continue
                    elif (not read.is_good(self.min_enddist, self.max_mismatch, self.mapq)) and (not mem_read.is_good(self.min_enddist, self.max_mismatch, self.mapq)):
                        #n_filtered += 1
                        n_filtered[X] += 1
                        continue
                    else:
                        n_nonfiltered[X] += 1

                    if X != R:
                        seen_alt.add(X)
                        n_pos[X] += 1
                        n_neg[X] += 1
        
                    #if not no_update:   
                    n_double[X] += 1

                    read_MQ = (read.mapq + mem_read.mapq)/2
                        
                    if read.base_qual > mem_read.base_qual:
                        read_BQ = (read.base_qual, mem_read.base_qual)
                        BQ_pair = f'({read.base_qual},{mem_read.base_qual})'
                    else:
                        read_BQ = (mem_read.base_qual, read.base_qual)
                        BQ_pair = f'({mem_read.base_qual},{read.base_qual})'

                    observed_BQs.add(BQ_pair)
                    enddist = max(read.enddist, mem_read.enddist)
                    has_indel = max(read.has_indel, mem_read.has_indel)
                    has_clip = max(read.has_clip, mem_read.has_clip)
                    NM = max(read.NM, mem_read.NM)
                    fragment_id = (min(read.start, mem_read.start), max(read.end, mem_read.end), X)
                    assert read.start < read.end
                    assert mem_read.start < mem_read.end
                    if fragment_id not in events[X]:
                        events[X][fragment_id] = (read_BQ, read_MQ, enddist, has_indel, has_clip, NM, BQ_pair)
                    else:
                        o_read_BQ = events[X][fragment_id][0]
                        if o_read_BQ < read_BQ:
                            events[X][fragment_id] = (read_BQ, read_MQ, enddist, has_indel, has_clip, NM, BQ_pair)

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
            if not read.is_good(self.min_enddist, self.max_mismatch, self.mapq):
                n_filtered[X] += 1
                continue
            else:    
                n_nonfiltered[X] += 1

            if X != R:
                seen_alt.add(X)
                if read.is_reverse:
                    n_neg[X] += 1
                else:
                    n_pos[X] += 1
            
            observed_BQs.add(str(read.base_qual))
            frag_id = (read.start, X)
            if frag_id not in events[X]:
                events[X][frag_id] = (read.base_qual, read.mapq, read.enddist, read.has_indel, read.has_clip, read.NM, str(read.base_qual))
            else:
                o_read_BQ = events[X][frag_id][0]
                if o_read_BQ < read.base_qual:
                    events[X][frag_id] = (read.base_qual, read.mapq, read.enddist, read.has_indel, read.has_clip, read.NM, str(read.base_qual))                        

        if len(seen_alt) == 0:
            return {}, {}, n_mismatch, n_double, n_pos, n_neg, n_filtered, n_nonfiltered


        new_mut_probs = defaultdict(dict)

        #I only need to calculate probabilities of changing bases from one of the seen alleles.
        # I have to considder change to all bases to calculate stay types (X->X) correctly.
        relevant_bases = [ref] + list(seen_alt)

        # Do positional update of error rates:
        for BQ in observed_BQs:
            new_mut_probs[BQ] = defaultdict(dict)

            for from_base in relevant_bases:
                p_rest = 1.0
                stay_type, stay_kmer = mut_type(from_base, from_base, ref_kmer)
                for to_base in ['A', 'C', 'G', 'T']:
                    if to_base == from_base:
                        continue
                    change_type, change_kmer = mut_type(from_base, to_base, ref_kmer)

                    p_prior = self.mut_probs[BQ][change_type][change_kmer]

                    a = p_prior * self.prior_N
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


        posterior_base_probs = {}
        BQs = {'A':[], 'C':[], 'G':[], 'T':[]}        


        for A in seen_alt:
            posterior_base_probs[A] = []
            for read_BQ, read_MQ, enddist, has_indel, has_clip, NM, str_BQ in events[A].values():
                muttype_from_A, kmer_from_A = mut_type(A, A, ref_kmer)
                muttype_from_R, kmer_from_R = mut_type(R, A, ref_kmer)
                posterior_from_A = new_mut_probs[str_BQ][muttype_from_A][kmer_from_A]
                posterior_from_R = new_mut_probs[str_BQ][muttype_from_R][kmer_from_R]

                posterior_base_probs[A].append((posterior_from_A, posterior_from_R, read_MQ))
                #if A==X:
                if type(read_BQ) == tuple:
                    BQ1, BQ2 = read_BQ
                    read_BQ = max(BQ1,BQ2)
                BQs[A].append((read_BQ, posterior_from_R, enddist, has_indel, has_clip, NM, str_BQ))

        for read_BQ, read_MQ, enddist, has_indel, has_clip, NM, str_BQ in events[R].values():
            muttype_from_A, kmer_from_A = mut_type(A, R, ref_kmer)
            muttype_from_R, kmer_from_R = mut_type(R, R, ref_kmer)
            posterior_from_A = new_mut_probs[str_BQ][muttype_from_A][kmer_from_A]
            posterior_from_R = new_mut_probs[str_BQ][muttype_from_R][kmer_from_R]
            for A in seen_alt:
                posterior_base_probs[A].append((posterior_from_A, posterior_from_R, read_MQ))
                #if A==X:
            if type(read_BQ) == tuple:
                BQ1, BQ2 = read_BQ
                read_BQ = max(BQ1,BQ2)
            BQs[R].append((read_BQ, posterior_from_R, enddist, has_indel, has_clip, NM, str_BQ))
        

        return posterior_base_probs, BQs, n_mismatch, n_double, n_pos, n_neg, n_filtered, n_nonfiltered


