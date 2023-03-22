# HERO
HERO is a hybrid error correction approach, which utilize short reads to correct long reads. HERO combines the merit of DBG and OG, but HERO addresses the shortage of OG by phasing reads with SNPs to avoid misalignment.

## Installation and dependencies
Please note that StrainXpress is built for linux-based systems and python3 only.
HERO relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)

To install HERO, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n hero
conda activate hero
conda install -c bioconda python=3.6 racon minimap2 fmlrc2 ratatosk
```
Subsequently, pull down the code to the directory and then you can directly use it:
```
git clone git@github.com:kangxiongbin/HERO.git

```
## Examples
We recommend first correcting the long reads with fmlrc2 or ratatosk.

Requires the input file to be uncompressed fastq or fasta. The short reads read1 and read2 must be in one file.
