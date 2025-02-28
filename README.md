# BBQ (BetterBaseQuals)

A somatic SNV variant caller that uses discordant overlapping read pairs to calculate sample-specific base error rates.

## Requirements

bbq requires Python 3.8 or above.

## Install

bbq can be installed using pip:
```
pip install betterbasequals
```
or using [pipx](https://pipx.pypa.io/stable/):
```
pipx install betterbasequals
```

### Usage

The bbq tool has different commands that can be run. To get at list of the different command you can use `-h` or `--help`:
```
bbq -h
```
To see the options for a specific command you can use `-h` after a command:
```
bbq count -h
```

To call variants in a sample you will need to first gather training data then train the model on the sample and finally call variants using the model. Typically you will want to split this into multiple jobs and we recommend using the snakemake workflow [bbq_pipeline](https://github.com/BesenbacherLab/bbq_pipeline) that does this. That github project also contain training data that can be used to test BBQ.

If you want to run both the error counting, model training and variant calling using a single command that is also possible using the following command (but we strongly recommend using the above mentioned workflow):

```
bbq call --bam_file cfdna.bam --filter_bam_file  pbmc.bam --twobit_file ~/Data/2bit/hg38.2bit
```

For some of the command it is possible to run it on only a part of the genome using the `-r {chrom}:{start}-{end}` option Fx:
```
bbq call_only --bam_file tmp.bam --twobit_file ~/Data/2bit/hg38.2bit --input_file_kmerpapa model.txt -r chr22:30000000-31000000 
```

Several of the commands requires the reference genome in 2bit format. If the bam file is mapped to a known reference genome the 2bit file can be downloaded from:
`https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.2bit` where `{genome}` is a valid UCSC genome assembly name (fx. "hg38").

Otherwise a 2bit file can be created from a fasta file using the `faToTwoBit` command. A linux binary of the tool can be downloaded from UCSC:
```
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
```

