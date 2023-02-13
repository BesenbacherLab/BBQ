import py2bit
from betterbasequals.utils import eprint, reverse_complement, Read, zip_pileups_single_chrom, open_bam_w_index, p2phred
from collections import defaultdict, Counter


def get_mut_type(ref, alt):
    if ref in ['T', 'G']:
        mtype = reverse_complement(ref) + '->' + reverse_complement(alt)
    else:
        mtype = ref + '->' + alt
    return mtype


class EmpiricalQuality:
    def __init__(
        self,
        bam_file,
        validation_bam_file,
        twobit_file,
    ):

        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)
        self.validation_bam_file = open_bam_w_index(validation_bam_file)
        self.tb = py2bit.open(twobit_file)
        self.radius = 2
        self.min_filter_depth = 10
        self.min_filter_BQ = 80

    def __del__(self):
        self.tb.close()

    def count_mutations_all_chroms(self):
        n_double = defaultdict(Counter)
        n_single = defaultdict(Counter)
        for idx_stats in self.bam_file.get_index_statistics():
            if idx_stats.mapped > 0:
                self.count_mutations(idx_stats.contig, None, None, n_double, n_single)
        return n_double, n_single


    def count_mutations(self, chrom, start, stop, n_double= None, n_single = None, mapq=50, mapq_hifi=40, prefix=""):
        pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            max_depth = 1000000,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=0,  # proper paired
            flag_filter=3840,
            min_base_quality = 1,
        )
        Hifi_pileup = self.validation_bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq_hifi,
            flag_filter=3840,
            min_base_quality = self.min_filter_BQ,
        )
        if n_double is None:
            n_double = defaultdict(Counter)
        if n_single is None:
            n_single = defaultdict(Counter)
        n_filtered_sites = 0
        n_used_sites = 0
        if start is None:
            eprint(f'Counting alleles at {chrom}.')
        else:
            eprint(f'Counting alleles at {chrom}:{start}-{stop}')

        for pileupcolumn, hifi_pc in zip_pileups_single_chrom(pileup, Hifi_pileup):
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
            
            hifi_alleles = set()
            N_hifi = 0
            for pread in hifi_pc.pileups:
                pos = pread.query_position
                if not pos is None:
                    hifi_alleles.add(pread.alignment.query_sequence[pos])
                    N_hifi += 1
            
            # We only considder good quality homo ref sites.
            if N_hifi < self.min_filter_depth or len(hifi_alleles) != 1:
                n_filtered_sites += 1
                continue
            allele = hifi_alleles.pop()
            if allele != ref:
                n_filtered_sites += 1
                continue

            reads_mem = {}

            coverage = 0

            muttype = {x: get_mut_type(ref, x) for x in ['A','C','G','T']}
            n_alt = 0
            for pileup_read in pileupcolumn.pileups:
                # test for deletion at pileup
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                coverage += 1
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
                events = []
                # Look for read partner
                if read.query_name in reads_mem:
                    # found partner process read pair
                    mem_read = reads_mem.pop(read.query_name)

                    if read.allel == mem_read.allel:
                        if read.allel != ref:
                            n_alt +=1
                        if read.base_qual != mem_read.base_qual:
                            continue
                        if read.isR1:
                            events.append(("double", read.base_qual, read.pos, True, muttype[read.allel]))
                            #n_double[(read.base_qual, read.pos, True)][muttype[read.allel]] += 1
                        else:
                            events.append(("double", read.base_qual, mem_read.pos, True, muttype[read.allel]))
                            #n_double[(read.base_qual, mem_read.pos, True)][muttype[read.allel]] += 1
                else:            
                    reads_mem[read.query_name] = read
            
            # Handle reads without partner (ie. no overlap)
            for read in reads_mem.values():
                if read.allel != ref:
                    n_alt += 1
                events.append(("single", read.base_qual, read.pos, read.isR1, muttype[read.allel]))
                #n_single[("single", read.base_qual, read.pos, read.isR1)][muttype[read.allel]] += 1


            if coverage > 10 and n_alt < 3:
                n_used_sites += 1
                for atype, BQ, pos, isR1, mtype in events:
                    if atype == "single":
                        n_single[("single", BQ, pos, isR1)][mtype] += 1
                    else:
                        n_double[("double", BQ, pos, isR1)][mtype] += 1
            else:
                n_filtered_sites += 1
        
        eprint(f'Filtered {n_filtered_sites} sites due to low coverage or alternative alleles in filter_bam_file.')
        eprint(f'Counted alleles at {n_used_sites} sites.')

        return n_double, n_single

