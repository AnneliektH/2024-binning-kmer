import csv
import argparse
from collections import defaultdict, OrderedDict

def read_manysearch(input_csv, msD, all_queries, contig_abunds):
    with open(input_csv, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            contig_name = row['match_name']
            query_name = row['query_name']
            all_queries.add(query_name)
            contig_abunds[contig_name].append(float(row['median_abund']))
            infoD = {f'{query_name}': float(row['median_abund']), f'{query_name}-var': float(row['std_abund']) **2}
            # this could get really slow with lots of samples. Could use tuple instead in that case...
            msD[contig_name].update(infoD)
    return msD, all_queries, contig_abunds

def write_results(output_tsv, allqueries, msData, contig_lengths, average_abunds):
    # build full list of fields (columns)
    fieldnames = ['contigName', 'contigLen', 'totalAvgDepth']
    query_columns = []
    for query_name in allqueries:
        query_columns.extend([f"{query_name}", f"{query_name}-var"])
    fieldnames.extend(query_columns)

    # initialize all with 0 to avoid missing cols getting set as NaN
    zeroes = {x: 0 for x in query_columns}

    with open(output_tsv, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        # results need to be in the same order as the input fasta, AND all contigs need to be represented, even if there was no abundance info
        for contig_name, contig_len in contig_lengths.items():
            # init all queries--> 0 to avoid NaNs
            contig_info = zeroes.copy()
            # update with contig info
            contig_info['contigName'] = contig_name
            contig_info['contigLen'] = contig_len
            total_avg_abund = average_abunds.get(contig_name, 0)
            contig_info['totalAvgDepth'] = total_avg_abund

            # get manysearch data for this contig, if we have any
            msInfoRow = msData.get(contig_name)
            if msInfoRow:
                contig_info.update(msInfoRow)
            writer.writerow(contig_info)


def main(args):

    msData = defaultdict(dict)
    allqueries = set()
    contig_abunds = defaultdict(list)

    # read the list of input CSV files
    input_files = [line.strip() for line in open(args.input_file_list, 'r')]

    # read contig lengths file
    contig_lengths = OrderedDict()
    with open(args.lengths, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            contig_lengths[row['name']] = int(row['length'])

    # process each file
    for input_manysearch in input_files:
        msData, allqueries, contig_abunds = read_manysearch(input_manysearch, msData, allqueries, contig_abunds)

    # get average abundance for each contig
    average_abunds = {}
    for contig_name, abunds in contig_abunds.items():
        average_abunds[contig_name] = sum(abunds) / len(abunds) if abunds else 0

    # Write final results
    write_results(args.output_tsv, allqueries, msData, contig_lengths, average_abunds)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate data from sourmash branchwater 'manysearch' to a metabat-formatted depth file.")
    parser.add_argument('input_file_list', type=str, help='File containing the list of input CSV files')
    parser.add_argument('-o', '--output_tsv', type=str, help='Path to the output TSV file', required=True)
    parser.add_argument('--lengths', type=str, help='Path to the output CSV file', required=True)

    args = parser.parse_args()
    main(args)
