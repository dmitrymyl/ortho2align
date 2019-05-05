import argparse
from subprocess import Popen
import pandas as pd
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Retrieve protein neighbourhood for given RNAs in bed format.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-rnafile',
                        type=str,
                        nargs='?',
                        help='file with RNAs in bed-like format.')
    parser.add_argument('-proteinfile',
                        type='str',
                        nargs='?',
                        help='protein annotation in gff/gtf/bed format.')
    parser.add_argument('-distance',
                        type=int,
                        nargs='?',
                        help='neighbourhood radius.')
    parser.add_argument('-rnaheaders',
                        type=str,
                        nargs='+',
                        help='headers of rna file.')
    parser.add_argument('-interfile',
                        type=str,
                        nargs='?',
                        default='protein_neighbourhood.tsv',
                        help='name of intermediate file.')
    parser.add_argument('-outfile',
                        type='str',
                        nargs='?',
                        default='protein_neighbourhood.json',
                        help='output filename.')
    parser.add_argument('-proteinid',
                        type=str,
                        nargs='?',
                        help='regex pattern to retrieve protein ids.')
    parser.add_argument('-rnaid',
                        type=str,
                        nargs='?',
                        help='column name with unique rna ids in rna file.')
    args = parser.parse_args()
    protein_headers = {"gff": ['prot_chr',
                               'prot_source',
                               'prot_cat',
                               'prot_start',
                               'prot_end',
                               'prot_smth1',
                               'prot_strand',
                               'prot_smth2',
                               'prot_data',
                               ],
                       "gtf": ['prot_chr',
                               'prot_source',
                               'prot_cat',
                               'prot_start',
                               'prot_end',
                               'prot_smth1',
                               'prot_strand',
                               'prot_smth2',
                               'prot_data',
                               ],
                       "bed": ['prot_chr',
                               'prot_start',
                               'prot_end',
                               'prot_data',
                               'prot_score',
                               'prot_strand',
                               ],
                       }

    find_neighbourhood = Popen(f"bedtools window -w {args.distance} -a {args.rnafile} -b {args.proteinfile} > {args.interfile}",
                               shell=True)
    find_neighbourhood.communicate()
    neighbourhood = pd.read_csv("protein_neighbourhood.tsv",
                                header=None,
                                index_col=None,
                                sep="\t")
    headers = args.rnaheaders + protein_headers[args.proteinfile.split(".")[-1]]
    neighbourhood.columns = headers
    prot_codes = pd.DataFrame(neighbourhood.loc[:, 'prot_data'].str.extract(args.proteinid))
    prot_codes.columns = ['prot_code']
    neighbourhood = pd.concat([neighbourhood, prot_codes], axis=1)
    rna_dicts = neighbourhood.loc[:, args.rnaheaders[0]:args.rnaheaders[-1]].drop_duplicates().to_dict('record')
    for i in range(len(rna_dicts)):
        current_rna = neighbourhood[neighbourhood.loc[:, args.rnaid] == rna_dicts[i][args.rnaid]]
        rna_dicts[i]['proteins'] = current_rna.loc[:, ['prot_chr',
                                                       'prot_start',
                                                       'prot_end',
                                                       'prot_strand',
                                                       'prot_code']
                                                   ].to_dict('record')
    with open(args.outfile, 'w') as outfile:
        json.dump(obj=rna_dicts, fp=outfile)
