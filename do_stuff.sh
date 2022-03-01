#!/bin/zsh
bamfile=/Volumes/GenomeDK/deep_sequencing_testis/DerivedData/bam/panTro5/Carl_Testis_150/Carl_Testis_150.aligned.sorted.markdup.bam 
filterbam=/Volumes/GenomeDK/MutationRates/NewVariantCalling/merged_bam_files/Carl.bam
twobit=~/Data/2bit/panTro5.2bit
poetry run betterbasequals  ${bamfile} ${filterbam} ${twobit} -r chr22:1-35000000 --output_file_good good_kmers.txt --output_file_bad bad_kmers.txt
