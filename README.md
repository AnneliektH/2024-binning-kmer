How to use sourmash to bin contigs on average abund of kmers, instead of using other bin algorithms?

HackMD: https://hackmd.io/muv86ruLR2aPq17e-72D0A?both

TEST DATA:
Atlas output of SRA dataset ERR3211994 (swine gut microbiome)
Atlas ran megahit, current min contig lenght is 1000 bp. Can lower to 200 if it increases binning.
ncontigs = 15,664

Atlas ran metabat and maxbin, and finalized binning using DAStool. 
Started out with 5 bins, based off the criteria left with 2:

Criteria (see folder quality_raport):
Completeness >50%
Contamination <5%
Length_scaffolds >=50000
Ambigious_bases <1e6
N50 > 5000
N_scaffolds < 1e3

Left with 2 bins: metabat_03 and metabat_09

TODO:
1. Calculate 3-mer frequencies (composition)
For each contig: count frequencies of 3-mers, merge reverse complements (n=32 3mers). Normalize by total number of 3mers in read

Sketch 3mers, singleton, scaled 50-100. 
build csv of frequencies of all 32 canonical 3-mers

2. Calculate coverage of 15-mers. *need to think about this more*
abundance sketch the 15 mers. 
sketch full file –> overall 15-mer abundances
stream through contigs/reads –> apply abundances from full file to k-mers found in contigs.

