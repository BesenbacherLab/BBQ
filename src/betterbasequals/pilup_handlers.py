from collections import Counter, defaultdict
from betterbasequals.utils import Read, ReadPair, reverse_complement, phred2p, p2phred

def get_mut_type(ref, papa_ref, alt):
    if ref != papa_ref:
        mtype = papa_ref + '->' + reverse_complement(alt)
    else:
        mtype = papa_ref + '->' + alt
    return mtype

ostrand = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
def mut_type(from_base, to_base, kmer):
    center = len(kmer)//2
    new_kmer = kmer[:center] + from_base + kmer[center+1:]
    if from_base in ['G','T']:
        return f'{ostrand[from_base]}->{ostrand[to_base]}', reverse_complement(new_kmer)
    else:
        return f'{from_base}->{to_base}', new_kmer


def get_pileup_count(pileupcolumn, ref, min_base_qual, blood=False):
    n_ref = 0
    n_alt = Counter()
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

        if not read.is_good():
            continue

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
    n_mismatch = Counter()
    n_double = 0

    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        #TODO: should consider what the right solution is if there is deletion at overlap

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

            n_double += 1
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
                adjusted_base_qual = max(adjusted_base_qual1, adjusted_base_qual2) + 3
                unadjusted_base_qual = max(read.base_qual, mem_read.base_qual)
                base_quals[read.allel].append((adjusted_base_qual, unadjusted_base_qual, 1))
            else:
                
                #Are ignoring ref alleles pt... should adjust later
                if read.allel != ref:
                    n_mismatch[read.allel] += 1
                    base_quals[read.allel].append((0, read.base_qual, 2))
                else:
                    n_ref += 1
                #Are ignoring ref alleles pt... should adjust later
                if mem_read.allel != ref:
                    n_mismatch[mem_read.allel] += 1
                    base_quals[mem_read.allel].append((0, mem_read.base_qual, 2))
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
    return base_quals, n_ref, n_mismatch, n_double


def get_alleles_w_probabities_seperate(pileupcolumn, ref, ref_kmer, correction_factor, improve = 1, ignore_ref = False, double_adjustment="max_plus3"):
    """
    Returns a dictionary that maps from allele A to a list of tuples with probability of 
    observing Alt allele A given then read and probability of observing ref allele R given 
    the read. I.e.: base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
    The probabilities are given on phred scale.
    We only considder reads where X_read_i == A or X_read_i == R.
    """

    reads_mem = {}
    base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    seen_alt = set()
    n_mismatch = Counter()
    n_double = Counter()
    n_alt = Counter()

    #tmp: counting ref... can be removed later
    n_ref = 0
    R = ref
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        #TODO: should consider what the right solution is if there is deletion at overlap

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

            # We do not trust mismathces in overlaps so we only handle matches
            if read.allel == mem_read.allel:
                X = read.allel
            
                if X == R:
                    if ignore_ref:
                        continue
                    alts = [A for A in ['A','C','G','T'] if A!=R]
                else:
                    alts = [X]
                    seen_alt.add(X)

                for A in alts:
                    n_double[A] += 1
                    
                    muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
                    a1, b1 = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
                    a2, b2 = correction_factor[mem_read.base_qual][muttype_from_A][kmer_from_A]
                    
                    if b1 is None:
                        adjusted_base_qual_from_A = ((a1+a2)/(2*improve), None)
                    else:
                        adjusted_base_qual_from_A = (a1+a2, improve*(b1+b2))
                    
                    muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
                    a1, b1 = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
                    a2, b2 = correction_factor[mem_read.base_qual][muttype_from_R][kmer_from_R]
                    
                    if b1 is None:
                        adjusted_base_qual_from_R = ((a1+a2)/(2*improve), None)
                    else:
                        #print("double", a1, a2, b1, b2, muttype_from_R, kmer_from_R)
                        adjusted_base_qual_from_R = (a1+a2, improve*(b1+b2))
                    BQ = max(read.base_qual, mem_read.base_qual)
                    if BQ > 35:
                        n_alt[X] += 1
                    base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R, BQ))

                    # adjusted_base_qual_from_A1 = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
                    # adjusted_base_qual_from_A2 = correction_factor[mem_read.base_qual][muttype_from_A][kmer_from_A]
                    
                    # if double_adjustment == "mult":
                    #     adjusted_base_qual_from_A = adjusted_base_qual_from_A1 + adjusted_base_qual_from_A2
                    # elif double_adjustment == "max_plus3":
                    #     adjusted_base_qual_from_A = max(adjusted_base_qual_from_A1, adjusted_base_qual_from_A2) + 3

                    # muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
                    # adjusted_base_qual_from_R1 = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
                    # adjusted_base_qual_from_R2 = correction_factor[mem_read.base_qual][muttype_from_R][kmer_from_R]
                    
                    # if double_adjustment == "mult":
                    #     adjusted_base_qual_from_R = adjusted_base_qual_from_R1 + adjusted_base_qual_from_R2
                    # elif double_adjustment == "max_plus3":
                    #     adjusted_base_qual_from_R = max(adjusted_base_qual_from_R1, adjusted_base_qual_from_R2) + 3
                        
                    # base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R))
            else: # Mismatch
                if read.allel != ref:
                    n_mismatch[read.allel] += 1
                    n_double[read.allel] += 1
                if mem_read.allel != ref:
                    n_mismatch[mem_read.allel] += 1
                    n_double[mem_read.allel] += 1
        else:            
            reads_mem[read.query_name] = read

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        X = read.allel
            
        if X == R:
            alts = [A for A in ['A','C','G','T'] if A!=R]
        else:
            alts = [X]
            seen_alt.add(X)
            if read.base_qual > 35:
                n_alt[X] += 1
        
        for A in alts:
            muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
            adjusted_base_qual_from_A = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
            
            muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
            adjusted_base_qual_from_R = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
            
            base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R, read.base_qual))
    
    posterior_base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    for A in base_probs:
        for (from_A, from_R, BQ) in base_probs[A]:
            a_from_A, b_from_A  = from_A
            a_from_R, b_from_R  = from_R 
            
            ## TODO: Should I do bayesian update on both from_A and from_R
            ## Or only on from_R ?????

            ## should I also save alpha and beta for posterior? 
            ## Or is it enough that I just save mean?

            if b_from_A is None:
                posterior_from_A = a_from_A
            else:
                posterior_from_A = (a_from_A + n_mismatch[A])/(a_from_A + b_from_A + n_double[A])

            if b_from_R is None:
                posterior_from_R = a_from_R
            else:
                posterior_from_R = (a_from_R + n_mismatch[A])/(a_from_R + b_from_R + n_double[A])

            posterior_base_probs[A].append((posterior_from_A, posterior_from_R))

    return base_probs, posterior_base_probs, seen_alt, n_mismatch, n_double, n_alt


# def update_BQ(BQ, N, n_mismatch, n_double):
#     p = phred2p(BQ)
#     a = p * N
#     b = N - a
#     new_p = (a + n_mismatch)/(a + b + n_double)
#     new_BQ = p2phred(new_p)
#     return new_BQ

# def get_alleles_w_probabities(pileupcolumn, ref, ref_kmer, correction_factor, N, double_adjustment="max_plus3"):
#     """
#     Returns a dictionary that maps from allele A to a list of tuples with probability of 
#     observing Alt allele A given then read and probability of observing ref allele R given 
#     the read. I.e.: base_probs[A] = [(P(A -> X_read_i|read_i),P(R -> X_read_i|read_i), ..., ]
#     The probabilities are given on phred scale.
#     We only considder reads where X_read_i == A or X_read_i == R.
#     """

#     reads_mem = {}
#     base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
#     seen_alt = set()
#     n_mismatch = Counter()
#     n_double = Counter()
#     n_alt = Counter()

#     #tmp: counting ref... can be removed later
#     n_ref = 0
#     R = ref
#     for pileup_read in pileupcolumn.pileups:
#         # test for deletion at pileup
#         if pileup_read.is_del or pileup_read.is_refskip:
#             continue
#         #TODO: should consider what the right solution is if there is deletion at overlap

#         # fetch read information
#         read = Read(pileup_read)

#         # test if read is okay
#         if (
#             read.allel not in "ATGC"
#             or read.start is None
#             or read.end is None
#             #or read.NH != 1
#         ):
#             continue


#         # Look for read partner
#         if read.query_name in reads_mem:
#             # found partner process read pair
#             mem_read = reads_mem.pop(read.query_name)

#             # We do not trust mismathces in overlaps so we only handle matches
#             if read.allel == mem_read.allel:
#                 X = read.allel
            
#                 if X == R:
#                     alts = [A for A in ['A','C','G','T'] if A!=R]
#                 else:
#                     alts = [X]
#                     seen_alt.add(X)

#                 for A in alts:
#                     n_double[A] += 1
                    
#                     muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
 
#                     adjusted_base_qual_from_A1 = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
#                     adjusted_base_qual_from_A2 = correction_factor[mem_read.base_qual][muttype_from_A][kmer_from_A]
                    
#                     if double_adjustment == "mult":
#                         adjusted_base_qual_from_A = adjusted_base_qual_from_A1 + adjusted_base_qual_from_A2
#                     elif double_adjustment == "max_plus3":
#                         adjusted_base_qual_from_A = max(adjusted_base_qual_from_A1, adjusted_base_qual_from_A2) + 3

#                     muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
#                     adjusted_base_qual_from_R1 = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
#                     adjusted_base_qual_from_R2 = correction_factor[mem_read.base_qual][muttype_from_R][kmer_from_R]
                    
#                     if double_adjustment == "mult":
#                         adjusted_base_qual_from_R = adjusted_base_qual_from_R1 + adjusted_base_qual_from_R2
#                     elif double_adjustment == "max_plus3":
#                         adjusted_base_qual_from_R = max(adjusted_base_qual_from_R1, adjusted_base_qual_from_R2) + 3
                        
#                     base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R))

#             else: # Mismatch
#                 if read.allel != ref:
#                     n_mismatch[read.allel] += 1
#                     n_double[read.allel] += 1
#                 if mem_read.allel != ref:
#                     n_mismatch[mem_read.allel] += 1
#                     n_double[mem_read.allel] += 1
#         else:            
#             reads_mem[read.query_name] = read

#     # Handle reads without partner (ie. no overlap)
#     for read in reads_mem.values():
#         X = read.allel
            
#         if X == R:
#             alts = [A for A in ['A','C','G','T'] if A!=R]
#         else:
#             alts = [X]
#             seen_alt.add(X)
#             if read.base_qual > 35:
#                 n_alt[X] += 1
        
#         for A in alts:
#             muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
#             adjusted_base_qual_from_A = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
            
#             muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
#             adjusted_base_qual_from_R = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
            
#             base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R))
    
#     posterior_base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
#     if N > 0:
#         for A in base_probs:
#             L = []
#             for (from_A, from_R) in base_probs[A]:
#                 L.append(update_BQ(from_A, N, n_mismatch[A], n_double[A]))

#             if b_from_A is None:
#                 posterior_from_A = a_from_A
#             else:
#                 posterior_from_A = (a_from_A + n_mismatch[A])/(a_from_A + b_from_A + n_double[A])

#             if b_from_R is None:
#                 posterior_from_R = a_from_R
#             else:
#                 posterior_from_R = (a_from_R + n_mismatch[A])/(a_from_R + b_from_R + n_double[A])

#             posterior_base_probs[A].append((posterior_from_A, posterior_from_R))

#     return base_probs, posterior_base_probs, seen_alt, n_mismatch, n_double, n_alt



def get_alleles_w_probabities_update(pileupcolumn, ref, ref_kmer, correction_factor, prior_N, no_update, double_adjustment=0.5, filter_reads = True):
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

                if filter_reads:
                    if (not read.is_good()) and mem_read.is_good():
                        #considder mem_read single read.
                        reads_mem[read.query_name] = mem_read
                        continue
                    elif (not mem_read.is_good()) and read.is_good():
                        #considder read single read.
                        reads_mem[read.query_name] = read
                        continue
                    elif (not read.is_good()) and (not mem_read.is_good()):
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
                    read_BQ = max(read.base_qual, mem_read.base_qual)
                    enddist = max(read.enddist, mem_read.enddist)
                    has_indel = max(read.has_indel, mem_read.has_indel)
                    has_clip = max(read.has_clip, mem_read.has_clip)
                    NM = max(read.NM, mem_read.NM)
                    #(read.base_qual, mem_read.base_qual)
                    events[A].append(("double", X, read_BQ, read_MQ, enddist, has_indel, has_clip, NM))

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
        if filter_reads and not read.is_good():
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
            events[A].append(("single", X, read.base_qual, read.mapq, read.enddist, read.has_indel, read.has_clip, read.NM))
        
    new_correction_factor = defaultdict(dict)

    #I only need to calculate probabilities of changing bases from one of the seen alleles.
    # I have to considder change to all bases to calculate stay types (X->X) correctly.
    relevant_bases = [ref] + list(seen_alt)

    for BQ in correction_factor:
        new_correction_factor[BQ]["single"] = defaultdict(dict)
        new_correction_factor[BQ]["double"] = defaultdict(dict)

        for from_base in relevant_bases:
            p_rest = 1.0
            p_rest_double = 1.0
            stay_type, stay_kmer = mut_type(from_base, from_base, ref_kmer)
            for to_base in ['A', 'C', 'G', 'T']:
                if to_base == from_base:
                    continue
                change_type, change_kmer = mut_type(from_base, to_base, ref_kmer)
                alpha, beta = correction_factor[BQ][change_type][change_kmer]
                p_prior = alpha/(alpha+beta)
                a = p_prior * prior_N
                b = prior_N - a
                if no_update or from_base != ref:
                    p_posterior = a/(a + b)
                else:
                    p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                p_rest -= p_posterior
                #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                new_correction_factor[BQ]["single"][change_type][change_kmer] = p2phred(p_posterior)

                p_prior_double = p_prior * double_adjustment
                a = p_prior_double * prior_N
                b = prior_N - a
                if no_update or from_base != ref:
                    p_posterior = a/(a + b)
                else:
                    p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                p_rest_double -= p_posterior
                new_correction_factor[BQ]["double"][change_type][change_kmer] = p2phred(p_posterior)
            
            new_correction_factor[BQ]["single"][stay_type][stay_kmer] = p2phred(p_rest)
            new_correction_factor[BQ]["double"][stay_type][stay_kmer] = p2phred(p_rest_double)

    posterior_base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    BQs = {'A':[], 'C':[], 'G':[], 'T':[]}
    for A in seen_alt:
        for overlap_type, X, read_BQ, read_MQ, enddist, has_indel, has_clip, NM in events[A]:
            muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
            muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
            posterior_from_A = new_correction_factor[read_BQ][overlap_type][muttype_from_A][kmer_from_A]
            posterior_from_R = new_correction_factor[read_BQ][overlap_type][muttype_from_R][kmer_from_R]

            # if overlap_type == "single":
            #     posterior_from_A = new_correction_factor[read_BQ][overlap_type][muttype_from_A][kmer_from_A]
            #     posterior_from_R = new_correction_factor[read_BQ][overlap_type][muttype_from_R][kmer_from_R]
            # elif overlap_type == "double":
            #     BQ1, BQ2 = read_BQ
            #     posterior_from_A = new_correction_factor[BQ1]["single"][muttype_from_A][kmer_from_A] + new_correction_factor[BQ2]["single"][muttype_from_A][kmer_from_A]
            #     posterior_from_R = new_correction_factor[BQ1]["single"][muttype_from_R][kmer_from_R] + new_correction_factor[BQ2]["single"][muttype_from_R][kmer_from_R]
            #     read_BQ = max(BQ1,BQ2)

            posterior_base_probs[A].append((posterior_from_A, posterior_from_R, read_MQ))
            if A==X:
                BQs[A].append((read_BQ, posterior_from_R, enddist, has_indel, has_clip, NM))
                    
    return posterior_base_probs, BQs, n_mismatch, n_double, n_pos, n_neg



def get_alleles_w_probabities_update_ver2(pileupcolumn, ref, ref_kmer, correction_factor, prior_N, no_update, double_adjustment=0.5):
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
    n_mismatch_BQ = {'A':Counter(), 'C':Counter(), 'G':Counter(), 'T':Counter()}
    n_double = Counter()
    n_double_BQ = {'A':Counter(), 'C':Counter(), 'G':Counter(), 'T':Counter()}
    n_alt = Counter()
    events = {'A':[], 'C':[], 'G':[], 'T':[]}

    R = ref
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        #TODO: should consider what the right solution is if there is deletion at overlap

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

            # We do not trust mismathces in overlaps so we only add to events list in case of match
            if read.allel == mem_read.allel:
                X = read.allel
            
                if X == R:
                    alts = [A for A in ['A','C','G','T'] if A!=R]
                else:
                    alts = [X]
                    seen_alt.add(X)

                for A in alts:
                    #if not no_update:    
                    n_double[A] += 1

                    read_MQ = (read.mapq + mem_read.mapq)/2
                    #if overlap_type == "double":
                    read_BQ = max(read.base_qual, mem_read.base_qual)
                    events[A].append(("double", X, read_BQ, read_MQ))

            else: # Mismatch
                #if not no_update:
                # TODO: Would it make sense to also count ref_mismatches so that we could
                # do bayesian update of A->R error rates and not only R->A error rates?
                if read.allel != ref and mem_read.allel == ref:
                    n_mismatch[read.allel] += 1
                    n_mismatch_BQ[read.allel][read.base_qual] += 1
                    for A in ['A','C','G','T']:
                        if A == ref:
                            continue
                        n_double[A] += 1
                        n_double_BQ[A][read.base_qual]
                if mem_read.allel != ref and read.allel == ref:
                    n_mismatch[mem_read.allel] += 1
                    n_mismatch_BQ[mem_read.allel][mem_read.base_qual] += 1
                    for A in ['A','C','G','T']:
                        if A == ref:
                            continue
                        n_double[A] += 1
                        n_double_BQ[A][mem_read.base_qual]
        else:            
            reads_mem[read.query_name] = read

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        X = read.allel
            
        if X == R:
            alts = [A for A in ['A','C','G','T'] if A!=R]
        else:
            alts = [X]
            seen_alt.add(X)
        
        for A in alts:
            events[A].append(("single", X, read.base_qual, read.mapq))
        
    new_correction_factor = defaultdict(dict)

    #I only need to calculate probabilities of changing bases from one of the seen alleles.
    # I have to considder change to all bases to calculate stay types (X->X) correctly.
    relevant_bases = [ref] + list(seen_alt)

    for BQ in correction_factor:
        new_correction_factor[BQ]["single"] = defaultdict(dict)
        new_correction_factor[BQ]["double"] = defaultdict(dict)

        for from_base in relevant_bases:
            p_rest = 1.0
            p_rest_double = 1.0
            stay_type, stay_kmer = mut_type(from_base, from_base, ref_kmer)
            for to_base in ['A', 'C', 'G', 'T']:
                if to_base == from_base:
                    continue
                change_type, change_kmer = mut_type(from_base, to_base, ref_kmer)
                alpha, beta = correction_factor[BQ][change_type][change_kmer]
                p_prior = alpha/(alpha+beta)
                a = p_prior * prior_N
                b = prior_N - a
                if no_update:
                    p_posterior = a/(a + b)
                else: 
                    p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                p_rest -= p_posterior
                #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                new_correction_factor[BQ]["single"][change_type][change_kmer] = p2phred(p_posterior)

                p_prior_double = p_prior * double_adjustment
                a = p_prior_double * prior_N
                b = prior_N - a
                if no_update:
                    p_posterior = a/(a + b)
                else:
                    p_posterior = (a + n_mismatch[to_base])/(a + b + n_double[to_base])
                #print(ref, BQ, from_base, to_base, p_posterior, n_mismatch[to_base], n_double[to_base])
                p_rest_double -= p_posterior
                new_correction_factor[BQ]["double"][change_type][change_kmer] = p2phred(p_posterior)
            
            new_correction_factor[BQ]["single"][stay_type][stay_kmer] = p2phred(p_rest)
            new_correction_factor[BQ]["double"][stay_type][stay_kmer] = p2phred(p_rest_double)

    posterior_base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    BQs = {'A':[], 'C':[], 'G':[], 'T':[]}
    for A in seen_alt:
        for overlap_type, X, read_BQ, read_MQ in events[A]:
            #muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
            muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
            
            #posterior_from_A = new_correction_factor[read_BQ][overlap_type][muttype_from_A][kmer_from_A]
            posterior_from_R = new_correction_factor[read_BQ][overlap_type][muttype_from_R][kmer_from_R]
            if A==X:
                posterior_base_probs[A].append((posterior_from_R, read_MQ, read_BQ))
            #if A==X:
            #    BQs[A].append((read_BQ, posterior_from_R))
                    
    return posterior_base_probs, n_mismatch, n_double, n_mismatch_BQ, n_double_BQ




def get_validation_probabities(pileupcolumn, ref, ref_kmer, correction_factor, improve = 1, ignore_ref = False, double_adjustment="max_plus3"):
    """
    
    """

    reads_mem = {}
    base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    seen_alt = set()
    n_mismatch = Counter()
    n_double = Counter()
    n_alt = Counter()

    #tmp: counting ref... can be removed later
    n_ref = 0
    R = ref
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        #TODO: should consider what the right solution is if there is deletion at overlap

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

            # We do not trust mismathces in overlaps so we only handle matches
            if read.allel == mem_read.allel:
                A = read.allel
                
                if A == R:
                    for X in [x for x in 'ACGT' if x != R]:
                        n_double[X] += 1
                    continue
                else:
                    n_double[A] += 1
                seen_alt.add(A)

                muttype_from_R, kmer_from_R = mut_type(R, A, ref_kmer)
                

                BQ = max(read.base_qual, mem_read.base_qual)

                alpha, beta = correction_factor[BQ][muttype_from_R][kmer_from_R]                
                
                if BQ > 35:
                    n_alt[A] += 1

                base_probs[A].append((alpha, beta, BQ, muttype_from_R, 1))

            else: # Mismatch
                if read.allel != ref:
                    n_mismatch[read.allel] += 1
                    n_double[read.allel] += 1
                if mem_read.allel != ref:
                    n_mismatch[mem_read.allel] += 1
                    n_double[mem_read.allel] += 1
        else:            
            reads_mem[read.query_name] = read

    # Handle reads without partner (ie. no overlap)
    for read in reads_mem.values():
        A = read.allel
                
        if A == R:
            continue

        seen_alt.add(A)
        if read.base_qual > 35:
            n_alt[A] += 1

        muttype_from_R, kmer_from_R = mut_type(R, A, ref_kmer)
        alpha, beta = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
        BQ = read.base_qual
        base_probs[A].append((alpha, beta, BQ, muttype_from_R, 2))
    
    #final_base_probs = {'A':[], 'C':[], 'G':[], 'T':[]}
    #for A in base_probs:
    #    for alpha, beta, BQ, muttype, atype in base_probs[A]:        
    #        final_base_probs[A].append((alpha, beta, BQ, muttype, atype, n_mismatch[A], n_))

    return base_probs, seen_alt, n_mismatch, n_double, n_alt





def get_adjustments(pileupcolumn, ref, papa_ref, kmer, correction_factor, change_dict):
    reads_mem = {}

    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        #TODO: should consider what the right solution is if there is deletion at overlap

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

                adjusted_base_qual1 = correction_factor[read.base_qual][mut_type][kmer]
                adjusted_base_qual2 = correction_factor[mem_read.base_qual][mut_type][kmer]

                #adjusted_base_qual1 = read.base_qual + correction_factor[mut_type][kmer]
                #adjusted_base_qual2 = mem_read.base_qual + correction_factor[mut_type][kmer]
                #adjusted_base_qual = adjusted_base_qual1 + adjusted_base_qual2

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
            adjusted_base_qual = correction_factor[read.base_qual][mut_type][kmer]
            change_dict[(read.query_name, read.isR1)].append((read.pos, read.base_qual, read.allel, int(adjusted_base_qual), 3))



def get_alleles_w_quals(pileupcolumn):
    base_quals = {'A':[], 'C':[], 'G':[], 'T':[]}
    for pileup_read in pileupcolumn.pileups:
        # test for deletion at pileup
        if pileup_read.is_del or pileup_read.is_refskip:
            continue

        # fetch read information
        read = Read(pileup_read)

        if not read.is_good():
            continue

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

