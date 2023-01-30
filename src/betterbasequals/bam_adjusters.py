import pysam
import py2bit
import sys
from betterbasequals.utils import reverse_complement, Read, open_bam_w_index
from betterbasequals.pilup_handlers import get_adjustments
from collections import defaultdict

class BaseAdjuster:
    def __init__(
        self,
        bam_file,
        twobit_file,
        kmer_papa,
        output_bam,
        adjustment_file,
        cutoff,
    ):
        # Open files for reading
        self.bam_file = open_bam_w_index(bam_file)

        #Open outputfile for writing
        self.out_file = pysam.AlignmentFile(output_bam, "wb", template=self.bam_file)

        self.tb = py2bit.open(twobit_file)

        #assumes that the kmer papa has been turned into phred scaled correction factor
        self.correction_factor = kmer_papa

        self.adjustment_file = adjustment_file
        self.cutoff = cutoff
                

    def __del__(self):
        self.tb.close()

    def call_mutations(self, chrom, start, stop, mapq=50, radius=3, prefix=""):
        pileup = self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
            flag_filter=3848,
        )

        change_dict = defaultdict(list)

        # Fill dict with changes:
        for pileupcolumn in pileup:
            ref_pos = pileupcolumn.reference_pos
            pileup_chrom = pileupcolumn.reference_name            

            try:
                kmer = self.tb.sequence(prefix + pileup_chrom, ref_pos- radius, ref_pos + radius + 1)
            except:
                continue
            
            ref = kmer[radius]
            if 'N' in kmer or len(kmer)!= 2*radius +1:
                continue            

            if kmer[radius] in ['T', 'G']:
                kmer = reverse_complement(kmer)
                papa_ref = reverse_complement(ref)
            else:
                papa_ref = ref

            assert len(ref) == 1
            if ref not in "ATGC":
                continue
            
            get_adjustments(pileupcolumn, ref, papa_ref, kmer, self.correction_factor, change_dict)

        print('n_reads_w_changes:', len(change_dict), file=sys.stderr)
        n_corrected_reads = 0
        n_uncorrected_reads = 0
        n_corrections = 0
        n_filtered = 0

        for read in self.bam_file.fetch(chrom, start, stop):
            if read.is_secondary or read.is_supplementary:
                n_filtered += 1
                continue
            read_id = (read.query_name, read.is_read1)
            if read_id in change_dict:
                n_corrected_reads += 1
                for pos, basequal, allele, adjusted_basequal, atype in change_dict[read_id]:
                    if not self.adjustment_file is None:
                        refpos = read.get_reference_positions(full_length=True)[pos]
                        print(read.reference_name, refpos, allele, basequal, adjusted_basequal, file = self.adjustment_file)
                    n_corrections += 1
                    assert(read.query_sequence[pos] == allele)
                    assert(read.query_qualities[pos] == basequal)
                    if self.cutoff is None or adjusted_basequal >= self.cutoff:
                        read.query_qualities[pos] = adjusted_basequal
                    else:
                        read.query_qualities[pos] = 0
            else:
                n_uncorrected_reads += 1
            self.out_file.write(read)
        return n_corrections, n_corrected_reads, n_uncorrected_reads, n_filtered

