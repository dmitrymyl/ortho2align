import json
import re
import argparse
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Transform given GFF/GTF annotation to json format.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-infile",
                        type=str,
                        nargs='?',
                        dest='infilename',
                        default='stdin',
                        help="input annotation filename, default stdin")
    parser.add_argument("-format",
                        type=str,
                        nargs='?',
                        dest='format',
                        choices=['gtf', 'gff'],
                        help="input annotation format")
    parser.add_argument("-outfile",
                        type=str,
                        nargs='?',
                        dest='outfilename',
                        help="output json filename")

    annotation_patterns = {'gff': ['chr',
                                   'source',
                                   'type',
                                   'start',
                                   'end',
                                   'smth1',
                                   'strand',
                                   'smth2',
                                   'data'],
                           'gtf': ['chr',
                                   'source',
                                   'type',
                                   'start',
                                   'end',
                                   'smth1',
                                   'strand',
                                   'smth2',
                                   'data']
                           }
    geneID_patterns = {'gff': r'GeneID:(\d+)',
                       'gtf': r'gene_id\s\"(\w*)\"'
                       }
    args = parser.parse_args()
    file_format = args.format
    if args.infilename == 'stdin':
        infile = sys.stdin
    else:
        infile = open(args.infilename, 'r')
    outfilename = args.outfilename
    annotation_pattern = annotation_patterns[file_format]
    geneID_pattern = geneID_patterns[file_format]
    holder = [dict(zip(annotation_pattern, line.strip().split("\t")))
              for line in infile
              if line[0] != "#"]
    if args.infilename != 'stdin':
        infile.close()
    keys = {re.search(geneID_pattern, entry['data']).group(1): entry
            for entry in holder
            if re.search(geneID_pattern, entry['data'])}
    with open(outfilename, 'w') as outfile:
        json.dump(keys, outfile)
