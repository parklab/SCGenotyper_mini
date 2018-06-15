# SCGenotyper_Mini

## Introduction
SCGenotyper_mini was written for the purpose of validating putative somatic mosaic SNVs using single cell DNA-seq data. It provides a lightweight solution to genotyping hundreds to thousands of sites across tens to hundreds of single-cell samples

For example, if you have a set of putative somatic SNVs discovered in deep bulk WGS and you would like to check within single-cell DNA-seq data for evidence of these events, then you're in the right place! SCGenotyper_mini uses local allele balance estimation to distinguish true variation from amplification artifact and sequencer error.

SCGenotyper_mini is not a new algorithm. Rather it is inspired by and draws heavily upon [SCcaller](https://github.com/biosinodx/SCcaller) with the following important differences:

  1. SCGenotyper_mini uses statistical haplotype-phase information to improve local allele balance estimation.
  2. SCGenotyper_mini uses a formal, parameterized likelihood-ratio test satisfying the condition of Wilk's Theorem to genotype positions rather than a bootstrapping method.
  3. SCGenotyper_mini **is not** a caller. While it will genotype a set of sites chosen a priori, it *will not* discover de novo sites of variation.

## Requirements
  1. Python 3.5 or greater
      * pandas
      * numpy
      * scipy
      * pysam (and necessarily samtools)
  
  2. A statistical haplotype phaser -- I recommend the [Eagle](https://data.broadinstitute.org/alkesgroup/Eagle/) software.
      * Alternatively, you can use a phasing / imputation server such as:
        * https://imputationserver.sph.umich.edu/
        * https://phasingserver.stats.ox.ac.uk/
    
## Usage

```
usage: run_mini.py [-h] [-b BAM [BAM ...]] [-B FBAM] -r POS -v VCF [-o FOUT]
                   [-l WINDOW] [--p-art P_ART] [--p-err P_ERR]
                   [--DP-min DP_MIN] [--mapq MAPQ]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM [BAM ...], --bam BAM [BAM ...]
                        bam(s) to genotype
  -B FBAM, --fbam FBAM  file containing paths to bam files to analyze.
                        Overrides -b.
  -r POS, --pos POS     tab-separated file of positions to genotype with
                        columns CHR, POS, REF, ALT
  -v VCF, --vcf VCF     vcf file of HAPLOTYPE-PHASED germline SNPs
  -o FOUT, --fout FOUT  output file (vcf format)
  -l WINDOW, --window WINDOW
                        window size for kernel
  --p-art P_ART         probability of amplification artifact
  --p-err P_ERR         probability of sequencer artifact
  --DP-min DP_MIN       Minimum depth to consider a training site
  --mapq MAPQ           Minimum mapq to consider a read
```

## Input Formats

### Inputing BAM files
You have two options for specifying BAM files to analyze. One must be specified or the program will raise an excpetion.
  1. `-b`: manually specify a list of bam files on the command line
  2. `-B`: a flat text file of BAM paths, one per line (supercedes `-b`)
  
### Inputing positions to genotype
Since SCGenotyper_mini *is not* a caller, you must provide a file specifing which positions to genotype. It should be a tap-separated file with columns `CHROM  POS REF ALT`. E.g.
```
12	16869088	G	T
10	28097311	A	G
1	27549935	A	C
1	196882197	A	G
7	112461481	G	T
```
The file should not include a header.

### Inputing germline SNPs
SCGentyper_mini requires a VCF of known germline SNP sites. Additionally, **this VCF must be haoplotype phased**. The easiest way to generate this file is to call SNPs using GATK's HaplotypeCaller from bulk sequencing of the same sample the single cells are derived from. Then phase the resulting VCF using Eagle either manually or by uploading to an imputation server listed above.
  * If phasing manually follow [Phasing with a refernce panel](https://data.broadinstitute.org/alkesgroup/Eagle/#x1-220005) of the Eagle manual.

## Output Formats
SCGentyper_mini produces a single, uncompressed VCF file with a row for every position specified with `-r` and a column for every sample specified with `-b` or `-B`.

The VCF uses the following genotypes:
```
./.: site cannot be genotyped
0|0: homozygous referenze
0|1: heterozygous with variant on the pseudo-maternal allele
1|0: heterozygous with a variant on the pseudo-paternal allele
1|1: homozygous alternate
0/.: homozygous reference with an amplification artifact
```

Each site in each sample is also assigned a filter flag:
```
NOSNP: no germline heterozyoug SNP within the kernel window (site cannot be genotyped)
LOWQUAL: The called genotype is not statistically distinguishable (p<0.05) from then next most likely genotype
ARTIF: The called genotype is not statistically distinguishable from an amplification artifact
PASS: All other filters are passed.
```
The FILTER flag for a site is set as follows: 
  * `PASS` if one or more samples has a `0|1` or `1|0` genotype which is `PASS` for that sample. 
  * `NOSNP` or `LOWQUAL` or `ARTIF` if all samples carry that respective flag
  * `.` otherwise.

Other format fields of interest are specified in the VCF header.

## FAQ
*1. A site is called 0|1 or 1|0 but no alternate reads are present. Why is this?*

When a position occurs within a region of allelic dropout (only one haplotype was sequenced), the estimated allele balance will be 0. Thus, in the absence of alternate reads, a heterozygous genotype is technically the most likely. However, since this genotype will be indistinguishable from `0|0` or potentially `0/.` it will always be flagged as `LOWQUAL` or `ARTIF`.
