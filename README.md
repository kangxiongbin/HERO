# HERO
HERO is a hybrid error correction approach that utilizes short reads to correct long reads. HERO combines the merits of the De Bruijn graph (DBG) and overlap graph (OG). However, it addresses the shortage of OG by phasing reads with SNPs to filter out the reads that originate from other strains/haplotypes and avoid overcorrection.

![workflow](https://github.com/kangxiongbin/HERO/blob/main/example/workflow.jpg=50%x)


The workflow of HERO. Black: long read r to be corrected. Yellow, blue, red: short reads, different colors indicate different haplotypes. Red crosses indicate errors. Tics (here: blue, red) indicate haplotype-specific variants.

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

Then HERO is used to further correct the pre-correct long reads with overlap-layout-consensus (OLC) paradigm. Requires the input file to be uncompressed fastq/fasta. The short reads read1 and read2 must be in one file. Additionally, there should not be spaces in the names of the reads, as this will affect the ability to correctly identify reads during downstream clustering. 

```
-Short reads correct long reads
python /folder/HERO/bin/HERO.py -r short_reads.fq -lc ratatosk3.fastq -p -o corrected_long.fa -s 30

-Long reads correct long reads
python /folder/HERO/bin/HERO.py -r ratatosk3.fastq -lc ratatosk3.fastq -o corrected_long.fa -i 1 -hlong_reads

```

More examples can be viewed and run on Code Ocean. 

```
https://codeocean.com/capsule/9666759/tree/v1
```

Note: If the reads files are quite large, for example Illumina reads files larger than 5G or long read files larger than 2G, it is recommended to split the files into smaller chunks for correction to improve speed. This can be achieved in HERO by setting -s to 200 or larger
