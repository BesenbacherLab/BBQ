from gwf import Workflow

gwf = Workflow(
    defaults={
        'cores': 1,
        'memory': '1g',
        'walltime': '11:00:00',
        'account' : 'MutationAnalysis'})


bamfile='/home/besen/mutationalscanning/DerivedData/bam/Illumina/panTro5/Carl_Testis_F003/Carl_Testis_F003.aligned.sorted.markdup.bam'
filterbam='/home/besen/MutationRates/NewVariantCalling/merged_bam_files/Carl.bam'
hifi_bam='/home/besen/mutationalscanning/DerivedData/bam/HiFi/panTro5/Carl_Testis_HiFi.bam'
twobit='/home/besen/Data/2bit/panTro5.2bit'
outfile = '/home/besen/mutationalscanning/Workspaces/besen/TestOverlap/results/BBQ_Carl_L3_1_chr22.txt.gz'
region='chr22'
gwf.target(
    f'test_chr22',
    inputs=[],
    outputs=[outfile]) << \
    f'''
    mkdir -p $(dirname {outfile})
    poetry run betterbasequals  {bamfile} {filterbam} {hifi_bam} {twobit} -r {region} --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt --output_file_kmerpapa kmer_papa_out.txt | gzip > {outfile}
    '''