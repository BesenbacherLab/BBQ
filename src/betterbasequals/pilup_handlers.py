from collections import Counter
from betterbasequals.utils import Read, ReadPair, reverse_complement

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


def get_alleles_w_probabities(pileupcolumn, ref, ref_kmer, correction_factor, double_adjustment="max_plus3"):
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
                    alts = [A for A in ['A','C','G','T'] if A!=R]
                else:
                    alts = [X]
                    seen_alt.add(X)

                for A in alts:
                    muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
                    adjusted_base_qual_from_A1 = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
                    adjusted_base_qual_from_A2 = correction_factor[mem_read.base_qual][muttype_from_A][kmer_from_A]
                    
                    if double_adjustment == "mult":
                        adjusted_base_qual_from_A = adjusted_base_qual_from_A1 + adjusted_base_qual_from_A2
                    elif double_adjustment == "max_plus3":
                        adjusted_base_qual_from_A = max(adjusted_base_qual_from_A1, adjusted_base_qual_from_A2) + 3

                    muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
                    adjusted_base_qual_from_R1 = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
                    adjusted_base_qual_from_R2 = correction_factor[mem_read.base_qual][muttype_from_R][kmer_from_R]
                    
                    if double_adjustment == "mult":
                        adjusted_base_qual_from_R = adjusted_base_qual_from_R1 + adjusted_base_qual_from_R2
                    elif double_adjustment == "max_plus3":
                        adjusted_base_qual_from_R = max(adjusted_base_qual_from_R1, adjusted_base_qual_from_R2) + 3
                        
                    base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R))

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
            muttype_from_A, kmer_from_A = mut_type(A, X, ref_kmer)
            adjusted_base_qual_from_A = correction_factor[read.base_qual][muttype_from_A][kmer_from_A]
            
            muttype_from_R, kmer_from_R = mut_type(R, X, ref_kmer)
            adjusted_base_qual_from_R = correction_factor[read.base_qual][muttype_from_R][kmer_from_R]
            
            base_probs[A].append((adjusted_base_qual_from_A, adjusted_base_qual_from_R))

    return base_probs, seen_alt



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

