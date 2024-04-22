import sys
import argparse
import csv
import screed

def main(args):

    with open(args.output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['name', 'length'])
        # parse through fasta records and write to csv
        records = screed.open(args.fastafile)
        for record in records:
            record_len = len(record.sequence)
            writer.writerow([record.name, record_len])

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("fastafile")
    p.add_argument("output_csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
