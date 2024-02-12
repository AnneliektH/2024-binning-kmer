import csv
import argparse
import screed
import sourmash

def main(args):
    """
    Report abundance of kmers in reads/contigs FROM whole file and output to csv
    """

    # first, load whole-file sketch and downsample if needed
    siglist = sourmash.load_file_as_signatures(args.input_sketch, ksize=args.ksize)
    siglist = list(siglist)
    assert len(siglist) == 1
    sig = siglist[0]
    fullfile_mh = sig.minhash.downsample(scaled=args.scaled)


    # generate fieldnames using hashvals
    fieldnames = ['name'] + list(fullfile_mh.hashes)

    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
    with open(args.output_csv, 'w', newline='') as csvfile:
        # Just use hashvals as column headers, shrug. 
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for record in screed.open(args.input_fasta):
            # first get existing hashvals
            mh.add_sequence(record.sequence)
            n_observed = len(mh) # total number kmers in this record
            # then inflate with fullsig counts but normalize by n kmers in the record
            kmers_to_ff_count = {hashval: (fullfile_mh.hashes[hashval]/n_observed) for hashval in mh.hashes}
            kmers_to_ff_count['name'] = record.name
            writer.writerow(kmers_to_ff_count)
            mh.copy_and_clear()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Get whole-file abundances for k-mers each record.')
    p.add_argument('input_fasta', type=str, help='Path to the input fasta file')
    p.add_argument('input_sketch', type=str, help='Path to corresponding sketch size')
    p.add_argument('-k', '--ksize', type=int, default=3, help='K-mer size (default: 3)')
    p.add_argument('-s', '--scaled', type=int, default=50, help='Scaled value for MinHash (default: 50)')
    p.add_argument('-o', '--output-csv', type=str, help='Path to the output CSV file')

    args = p.parse_args()
    main(args)
