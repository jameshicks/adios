# adios
ADIOS: Adaptive Detection of IBD over Sequence data

adios is a hidden-markov method for efficiently detecting genome segments shared IBD between pair of individuals from whole genome sequencing data.
The majority of genotypes between two individuals are not informative of IBD status. 
adios considers only the informative genotype pairs between two individuals: opposite homozygotes (excludes IBD) or loci where at least one individual has a rare variant.
Segments are returned with a LOD-score style quality metric to quantify the confidence in the IBD call.

## Options
+ `--vcf`: VCF input file
+ `--vcf_freq`: INFO field in vcf to use as allele frequency
+ `--keep_singletons`: Include singleton variants in dataset
+ `--keep_monomorphic`: Include monomorphic positions in dataset
+ `--empirical_freqs`: Calculate allele frequencies from data (supercedes `--vcf_freq`)
+ `--freq_floor`: Minimum frequency used in calculations
+ `--rare`: Alleles below this frequency are included as rare alleles.
+ `--minlength`: Minimum segment length to report (in megabases)
+ `--minmark`: Minimum number of markers to declare IBD
+ `--minlod`: Minimum segment LOD to declare IBD
+ `--err`: Genotype error rate.
+ `--transition`: Transition cost to enter and leave IBD states. (Two integers required, larger numbers correspond to higher cost)


## Installation 
adios can be generated by running `make`. Building adios currently requires a compiler that supports the C++11 standard (g++ 4.8 or clang 3.3). 
No external libraries are required. Running `make test` will run unit tests on many of the functions used (requires the library cpputest).

## Output:

1. IND_1: The first individual in the IBD pair 
2. IND_2: The first individual in the IBD pair
3. CHROM: Chromome
4. START: Segment start (bp)
5. END: Segment stop (bp)
6. LENGTH: Segment length (with unit)
7. STATE: IBD state (1 or 2)
8. NMARK: Number of markers used to declare IBD
9. NERR: Number of opposite homozygotes in segment
9. LOD: log10(L(segment|IBD=STATE) / L(segment|IBD=0)) 
