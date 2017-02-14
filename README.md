# adios
ADIOS: Adaptive Detection of IBD over Sequence data

adios is a hidden-markov method for efficiently detecting genome segments shared IBD between pair of individuals from whole genome sequencing data.
The majority of genotypes between two individuals are not informative of IBD status. 
adios considers only the informative genotype pairs between two individuals: opposite homozygotes (excludes IBD) or loci where at least one individual has a rare variant.
Segments are returned with a LOD-score style quality metric to quantify the confidence in the IBD call.

## Installation 
adios can be installed from source in the standard way: `./configure && make`. 
The configure script will detect the availability of OpenMP and, if present, enable multithreading.  
Building adios currently requires a compiler that supports the C++11 standard (g++ >= 4.8 or >= clang 3.3). 
No external libraries are required. Running `make test` will run unit tests on many of the functions used (requires the library cpputest).

## Options
+ `--vcf`: VCF input file
+ `--vcf_freq`: INFO field in vcf to use as allele frequency, otherwise calculated from data
+ `--out`: Prefix for output file
+ `--keep_singletons`: Include singleton variants in dataset
+ `--keep_monomorphic`: Include monomorphic positions in dataset
+ `--freq_floor`: Minimum frequency used in calculations
+ `--rare`: Alleles below this frequency are included as rare alleles.
+ `--minlength`: Minimum segment length to report (in megabases)
+ `--minmark`: Minimum number of shared rare variants to declare IBD
+ `--minlod`: Minimum segment LOD to declare IBD
+ `--err`: Genotype error rate.
+ `--transition`: Transition cost to enter and leave IBD states. (Two integers required, larger numbers correspond to higher cost)
+ `--threads`: Number of threads to split analysis over.


## Output:

1. IND_1: The first individual in the IBD pair 
2. IND_2: The first individual in the IBD pair
3. CHROM: Chromome
4. START: Segment start (bp)
5. END: Segment stop (bp)
6. LENGTH: Segment length (with unit)
7. STATE: IBD state (1 or 2)
8. NMARK: Number of informative markers in IBD segment
9. NRARE: Number of rare variants shared by both individuals 
10. NERR: Number of opposite homozygotes in segment
11. LOD: log10(L(segment|IBD=STATE) / L(segment|IBD=0)) 
