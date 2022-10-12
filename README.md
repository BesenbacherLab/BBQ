# BetterBaseQuals

Calculates sample-specific base qualities using overlapping reads.

## Requirements

BetterBaseQuals requires Python 3.7 or above.

## Development:

For development we use [poetry](https://python-poetry.org) to handle dependencies.
If poetry is not installed on your system you can install poetry using this command:
```
curl -sSL https://install.python-poetry.org | python3 -
```
Once poetry is installed you can install the current project including all dependencies using the following command:
```
poetry install
```

### Using poetry to run the software

poetry creates a virtual environment so that the current project is isolated from the rest of the system. To run a command in the poetry virtual environment you type `poetry run {command}`. So to run the bbq tool you can type:
```
poetry run bbq
```
This only works if you are in the current directory. If you want to run the code from a different directory you can first make sure that the current version is installed in the poetry environment by running `poetry install` and then using the `which` command to get the full path of the installed version:
```
poetry install
poetry run which bbq
```
That will return a full path to an executable that can be run from any directory.


### Using poetry to add dependencies

If we want to use numpy in our program we can add it as a dependency using the command:
```
poetry add numpy
```
poetry will then find a version of numpy that fits with the requirements of other packages and install it.


### Usage

The bbq tool has different commands that can be run. To get at list of the different command you can use `-h` or `--help`:
```
poetry run bbq -h
```
To see the options for a specific command you can use `-h` after a command:
```
poetry run bbq count -h
```

For some of the command it is possible to run it on only a part of the genome using the `-r {chrom}:{start}-{end}` option Fx:
```
poetry run adjust_only --bam_file sample.bam --twobit_file hg38.2bit --input_file_kmerpapa model.txt --outbam sample_adjusted.bam -r chr22:30000000-31000000 
```

Several of the commands requires the reference genome in 2bit format. If the bam file is mapped to a known reference genome the 2bit file can be downloaded from:
`https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.2bit` where `{genome}` is a valid UCSC genome assembly name (fx. "hg38").

Otherwise a 2bit file can be created from a fasta file using the `faToTwoBit` command. A linux binary can of the tool can be downloaded from UCSC:
```
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
```

