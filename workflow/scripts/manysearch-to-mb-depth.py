import csv
import argparse
from collections import defaultdict

def read_manysearch(input_csv, msD, all_queries, contig_abunds):
    with open(input_csv, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            contig_name = row['match_name']
            query_name = row['query_name']
            all_queries.add(query_name)
            contig_abunds[contig_name].append(float(row['average_abund']))
            infoD = {f'{query_name}': float(row['average_abund']), f'{query_name}-var': float(row['std_abund']) **2}
            # this could get really slow with lots of samples. Could use tuple instead in that case...
            msD[contig_name].update(infoD)
    return msD, all_queries, contig_abunds

def write_results(output_csv, allqueries, msData):
    # build full list of fields (columns)
    fieldnames = ['contigName', 'contigLen', 'totalAvgDepth']
    query_columns = []
    for query_name in allqueries:
        query_columns.extend([f"{query_name}", f"{query_name}-var"])
    fieldnames.extend(query_columns)
    
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

        for msInfo in msData.values():
            writer.writerow(msInfo)


def main(args):
    
    msData = defaultdict(dict)
    allqueries = set()
    contig_abunds = defaultdict(list)

    # read the list of input CSV files
    input_files = [line.strip() for line in open(args.input_file_list, 'r')]

    # read contig lengths file 
    contig_lengths = {}
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
    
    # store the average abundance for each contig in msData
    for contig_name, avg_abund in average_abunds.items():
        # fill out remaining columns for each contig
        msData[contig_name]['totalAvgDepth'] = avg_abund
        msData[contig_name]['contigName'] = contig_name
        if contig_name in contig_lengths:
            msData[contig_name]['contigLen'] = contig_lengths[contig_name]
        else:
            msData[contig_name]['contigLen'] = 0
            print(f"WARNING: contig {contig_name} not found in contig lengths file.")
         
    # Write final results
    write_results(args.output_csv, allqueries, msData)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate data from sourmash branchwater 'manysearch' to a metabat-formatted depth file.")
    parser.add_argument('input_file_list', type=str, help='File containing the list of input CSV files')
    parser.add_argument('-o', '--output_csv', type=str, help='Path to the output CSV file', required=True)
    parser.add_argument('--lengths', type=str, help='Path to the output CSV file', required=True)

    args = parser.parse_args()
    main(args)