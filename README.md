## Binning using sourmash mgmanysearch output
- Both VAMB and metabat2 use the jgi contig depth file for binning, which is created from read mapping 
- input for jgi_depth is bam files. 


Can we use sourmash output instead, where we use kmer abundances in each contig instead of mapping depths? Because read mapping takes forever.

