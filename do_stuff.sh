#!/bin/zsh
bamfile=/Volumes/GenomeDK/mutationalscanning/DerivedData/bam/Illumina/panTro5/Carl_Testis_F003/Carl_Testis_F003.aligned.sorted.markdup.bam
filterbam=/Volumes/GenomeDK/MutationRates/NewVariantCalling/merged_bam_files/Carl.bam
hifi_bam=/Volumes/GenomeDK/mutationalscanning/DerivedData/bam/HiFi/panTro5/Carl_Testis_HiFi.bam
twobit=~/Data/2bit/panTro5.2bit

#region=chr22:0-35000000
region=chr22:30000000-31000000
poetry run betterbasequals  ${bamfile} ${filterbam} ${hifi_bam} ${twobit} -r ${region} --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt --output_file_kmerpapa kmer_papa_out.txt 
#poetry run betterbasequals  ${bamfile} ${filterbam} ${twobit} -r chr22:30000000-30100000 --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt
