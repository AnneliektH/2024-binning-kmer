# srun
srun --account=ctbrowngrp -p bmm -J checm2 -t 4:00:00 -c 16 --mem=50gb --pty bash

# snakemake
mamba activate snakemake
snakemake --use-conda --resources mem_mb=25000 --rerun-triggers mtime -c 16 --rerun-incomplete -k 

# seqtk for cutting sequences from multifasta
mamba activate seqtk
seqtk subseq test.fa test.txt

for f in *.txt
do
seqtk subseq ../atlas_ERR3211994/ERR3211994/ERR3211994_contigs.fasta $f > $f.fa
ddone

# run checkm2 to check the 'bins'
srun --account=ctbrowngrp -p bmm -J checm2 -t 4:00:00 -c 16 --mem=50gb --pty bash
mamba activate checkm2
checkm2 predict --input ./fa_files/ --output-directory checkm_out -x fa --threads 16

# sketch contigs and readsfile 
# Contigs no abundance, read file with abundance
mamba activate sourmash 
sourmash sketch dna --singleton \
ERR3211994_contigs.fasta -p k=15,scaled=100  \
-o ../../sourmash_sig/ERR3211994.contigs.single.k15.sig.gz

sourmash sketch dna ERR3211994_clean_R1.fastq.gz -p abund,k=15,scaled=100 -o ../../../sourmash_sig/ERR3211994.R1.abund.k15.reads.sig.gz
sourmash sketch dna ERR3211994_clean_R2.fastq.gz -p abund,k=15,scaled=100 -o ../../../sourmash_sig/ERR3211994.R2.abund.k15.reads.sig.gz

# make the metabath depth file so we see what it needs to look like
mamba activate metabat
jgi_summarize_bam_contig_depths --percentIdentity 97 --outputDepth ../depth.metabat.txt ./ERR3211994/sequence_alignment/ERR3211994.bam 

A file having mean and variance of base coverage depth (tab delimited; the first column should be contig names, and the first row will be considered as the header and be skipped)
So contiglength, mean depth, depth variance.

# sourmash mgsearch, one metag per file, only flat sigs 
# change the source code to do the ksize and scale i want
in /home/amhorst/mambaforge/envs/sourmash/lib/python3.10/site-packages/sourmash_plugin_containment_search.py
sourmash scripts mgsearch ERR3211994_0.k15.sig ERR3211994.abund.k15.reads.sig.gz \
-k 15 --scaled 100

sourmash scripts mgmanysearch --queries ERR3211994.contigs.single.k15.sig.gz \
--against ERR3211994.R1.abund.k15.reads.sig.gz -k 15 --scaled 100 -o mgmanysearch.abund.csv


# try for one bin
sourmash scripts mgmanysearch --queries bin_3.sing.k15.sig.gz --against ERR3211994.R1.abund.k15.reads.sig.gz -k 15 --scaled 100 -o bin3.abund.csv

# cut bin into tiny pieces and try with that
