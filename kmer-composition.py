import csv
import argparse
import screed
from itertools import product
import sourmash
from sourmash.minhash import hash_murmur
import pytest


def canonical_kmers_and_hashes(ksize):
    """
    Generate all unique canonical k-mers for a specified k-mer size and their hash values.
    
    Parameters:
    - ksize: The size of the k-mer (int).
    
    Returns:
    A dictionary where keys are canonical k-mer strings and values are their hash values for the given k-mer size.
    """
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=ksize)]
    print(f'kmers at ksize {ksize}: ', len(all_kmers))
    hashes_to_kmers = {}
    seen_kmers=set()
    for kmer in all_kmers:
        # first, get canonical kmer
        rev_comp = screed.rc(kmer) # get reverse complement
        canonical_kmer = min(kmer, rev_comp)
        
        # Ensure we only compute and store the hash for each unique canonical k-mer once
        if canonical_kmer not in seen_kmers:
            hashval = hash_murmur(canonical_kmer)
            hashes_to_kmers[hashval] = canonical_kmer
            seen_kmers.add(canonical_kmer)
    print(f'canonical kmers at ksize {ksize}: ', len(seen_kmers))
    return hashes_to_kmers 


def main(args):
    """
    Count abundance of kmers in reads/contigs and output to csv
    """
    hashes_to_kmers = canonical_kmers_and_hashes(args.ksize)
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled, track_abundance=True)
    with open(args.output_csv, 'w', newline='') as csvfile:
        # Use kmer idents rather than hashvals as fieldnames -- nicer csv headers :)
        fieldnames = ['name'] + list(hashes_to_kmers.values())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for record in screed.open(args.input_fasta):
            # use sourmash to count abundances
            mh.add_sequence(record.sequence)
            # map hashvals to kmer idents to allow nicer csv headers :)
            n_observed = len(mh.hashes) # total number kmers in this record
            print(n_observed) # gives 0
            if n_observed:
                kmer_to_count = {hashes_to_kmers[hashval]: (mh.hashes.get(hashval, 0)/n_observed) for hashval in hashes_to_kmers}
                kmer_to_count['name'] = record.name
                writer.writerow(kmer_to_count)
            mh.copy_and_clear()

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Get kmer composition of each sequence record in a file')
    p.add_argument('input_fasta', type=str, help='Path to the input fasta file')
    p.add_argument('-k', '--ksize', type=int, default=3, help='K-mer size (default: 3)')
    p.add_argument('-s', '--scaled', type=int, default=50, help='Scaled value for MinHash (default: 50)')
    p.add_argument('-o', '--output-csv', type=str, help='Path to the output CSV file')

    args = p.parse_args()
    main(args)


# test canonical kmers fn
def test_canonical_kmers_and_hashes():
    kh = canonical_kmers_and_hashes(ksize=3)
    assert len(kh) == 32
    kh = canonical_kmers_and_hashes(ksize=4)
    assert len(kh) == 136
