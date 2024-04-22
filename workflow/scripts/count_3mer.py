# imports
import sourmash
import sys
import csv
from collections import Counter
import screed
import pandas as pd

# calculate reverse complement of a kmer, and return the canonical one (alphabetically smaller)
def canonkmer(kmer):
    rc_kmer = screed.rc(kmer)
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
    return canonical_kmer  


# count canonical kmers in a sequence, return counts
def count3mers(seq):
    counts=Counter()
    mh = sourmash.MinHash(n=0, ksize=3, scaled=10)
    for kmer, hashval in mh.kmers_and_hashes(seq):
        counts[(canonkmer(kmer))] += 1
    return counts


# Set the headers for csv file. For longer kmers it needs to be changed. 
all_3mers = [f'{a}{b}{c}' for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']
can_3mers = []
for mer in all_3mers:
    mer = canonkmer(mer)
    can_3mers.append(mer)
# keep only canonical headers
can_3mers = list(set(can_3mers))


# write counted kmers from each fasta to csv
def main(input_fasta, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['read_name'] + can_3mers
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
            
        for record in screed.open(input_fasta):
            kmer_counts = count3mers(record.sequence)
            row = {'read_name': record.name}
            row.update(kmer_counts)
            writer.writerow(row)
    df = pd.read_csv(output_csv)
    df = df.set_index(df['read_name'])
    del df['read_name']
    row_sums = df.sum(axis=1)
    normalized_df = df.div(row_sums, axis=0)
    normalized_df.to_csv(output_csv)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 3mercounter.py input.fasta output.csv")
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]
    main(input_fasta, output_csv)





