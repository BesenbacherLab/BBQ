from collections import Counter
from betterbasequals.utils import Read, ReadPair


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
