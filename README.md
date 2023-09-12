# HERO
HERO is a hybrid error correction approach that utilizes short reads to correct long reads or polish contigs. HERO combines the merits of the De Bruijn graph (DBG) and overlap graph (OG). However, it addresses the shortage of OG by phasing reads with SNPs to filter out the reads that originate from other strains/haplotypes and avoid overcorrection.

![overcorrection](https://github.com/kangxiongbin/HERO/assets/23208764/8d334ac1-e0f6-4ce3-b49a-381c70fad6f4)

As shown in the above figure, contigs represent unpolished contigs. The previous approach with Racon was to align short reads to contigs using minimap2 and then generate new consensus contigs with poa to achieve error correction. Similar to Racon, HERO also aligns reads to contigs using minimap2, but HERO utilizes SNP information to filter out short reads originating from other strains/haplotypes first, before calling Racon's POA to generate the consensus contigs. This avoids interference from incorrectly aligned reads. For example at position A, HERO can make the correction, but it will not introduce errors like BCDEFG as Racon does, which increases the mismatch error rates.

## Installation and dependencies
Please note that HERO is built for linux-based systems and python3 only.
HERO relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)

To install HERO, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n hero
conda activate hero
conda install -c bioconda python=3.6 racon minimap2 ropebwt2 fmlrc2 ratatosk lordec
```
Subsequently, pull down the code to the directory and then you can directly use it:
```
git clone git@github.com:kangxiongbin/HERO.git
```
## Examples
We recommend first correcting the long reads 3 times with DBG approaches such as ratatosk, fmlrc2 and lordec.

```
Ratatosk correct -v -c 30  -s short_reads.fq -l long_reads.fq -o ratatosk1
Ratatosk correct -v -c 30  -s short_reads.fq -l ratatosk1.fastq -o ratatosk2
Ratatosk correct -v -c 30  -s short_reads.fq -l ratatosk2.fastq -o ratatosk3
```

Then HERO is used to further correct the pre-correct long reads with overlap-layout-consensus (OLC) paradigm. Requires the input file to be uncompressed fastq/fasta. The short reads read1 and read2 must be in one file.

```
-Short reads correct long reads
python /folder/HERO/bin/HERO.py -r short_reads.fq -lc ratatosk3.fastq -p -o corrected_long.fa -s 30

-Long reads correct long reads
python /folder/HERO/bin/HERO.py -r ratatosk3.fastq -lc ratatosk3.fastq -o corrected_long.fa -i 1 -hlong_reads

```
