#!/bin/zsh
bamfile=/Volumes/GenomeDK/mutationalscanning/DerivedData/bam/Illumina/panTro5/Carl_Testis_F003/Carl_Testis_F003.aligned.sorted.markdup.bam
filterbam=/Volumes/GenomeDK/MutationRates/NewVariantCalling/merged_bam_files/Carl.bam
hifi_bam=/Volumes/GenomeDK/mutationalscanning/DerivedData/bam/HiFi/panTro5/Carl_Testis_HiFi.bam
#twobit=~/Data/2bit/panTro5.2bit
twobit=/Volumes/GenomeDK/Data/2bit/panTro5.2bit

#region=chr22:0-35000000
region=chr22:30000000-30100000
#poetry run betterbasequals  --bam${bamfile} ${filterbam} ${hifi_bam} ${twobit} -r ${region} --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt --output_file_kmerpapa kmer_papa_out.txt 
#poetry run betterbasequals  ${bamfile} ${filterbam} ${twobit} -r chr22:30000000-30100000 --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt
#poetry run betterbasequals 
echo /Users/au149884/Library/Caches/pypoetry/virtualenvs/betterbasequals-bpN7EUlj-py3.9/bin/betterbasequals --bam_file ${bamfile} --filter_bam_file ${filterbam} --validation_bam_file ${hifi_bam} --twobit_file ${twobit} -r ${region} --output_file_good good_kmers_chr10.txt --output_file_bad bad_kmers_chr10.txt --output_file_kmerpapa kmer_papa_out_chr10.txt

#/Users/au149884/Library/Caches/pypoetry/virtualenvs/betterbasequals-bpN7EUlj-py3.9/bin/betterbasequals --bam_file ${bamfile} --filter_bam_file ${filterbam} --validation_bam_file ${hifi_bam} --twobit_file ${twobit} -r ${region} --input_file_good good_kmers_chr10.txt --input_file_bad bad_kmers_chr10.txt --output_file_kmerpapa kmer_papa_out_chr10.txt --train_kmerpapas_only  

/Users/au149884/Library/Caches/pypoetry/virtualenvs/betterbasequals-bpN7EUlj-py3.9/bin/bbq count --bam_file /Volumes/GenomeDK/mutationalscanning/DerivedData/bam/Illumina/panTro5/Carl_Testis_F003/Carl_Testis_F003.aligned.sorted.markdup.bam --twobit_file /Volumes/GenomeDK/Data/2bit/panTro5.2bit -r chr22:34000000-35000000 --output_file_good good_kmers_chr10.txt --output_file_bad bad_kmers_chr10.txt

/Users/au149884/Library/Caches/pypoetry/virtualenvs/betterbasequals-bpN7EUlj-py3.9/bin/bbq train_only  --input_file_good good_kmers_chr10.txt --input_file_bad bad_kmers_chr10.txt --output_file_kmerpapa kmerpapa_chr22-34-35.txt
/Users/au149884/Library/Caches/pypoetry/virtualenvs/betterbasequals-bpN7EUlj-py3.9/bin/bbq test_kmerpapa  --input_file_good good_kmers_chr10.txt --input_file_bad bad_kmers_chr10.txt --input_file_kmerpapa kmerpapa_chr22-34-35.txt
