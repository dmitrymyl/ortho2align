import json
import os
import argparse
from subprocess import Popen, PIPE


dictfilt = lambda x, y: {k: v for k, v in x.items() if k in y}
get_coord = lambda x, y: "{}:{}-{}({})".format(x[y[0]], x[y[1]], x[y[2]], x.get(y[3], "."))


def get_bed(list_of_dicts, chrom='chr', start='start', end='end', strand='strand'):
    """
    Return bed-formated string for dictionary records of
    genomic ranges.

    Args:
        list_of_dicts: contatins genomic ranges coordinates.
        chrom (str): key name of chromosome in list_of_dicts records
            (default: 'chr').
        start (str): key name of start in list_of_dicts records
            (default: 'start').
        end (str): key name of end in list_of_dicts records
            (default: 'end').
        strand (str): key name of strand in list_of_dicts records
            (default: 'strand').

    Returns:
        bed-formated string with one genomic range per line.
    """
    data = ("{}\t{}\t{}\t.\t.\t{}".format(record[chrom], record[start],
                                          record[end], record.get(strand, "."))
            for record in list_of_dicts)
    return "\n".join(set(data))


def get_fasta(genome, bed_ranges, fname, directory):
    """
    Retrives sequences in fasta files from given genome
    corresponding to given bed ranges.

    Args:
        genome (str): path to corresponding genome.
        bed_ranges (str): bed-formatted string with
            one genomic range per line.
        fname (str): the name of intermediate dump file
            with all fasta sequences.
        directory (str): a path where to put extracted
            sequences.

    Returns:
        None. Each sequence corresponding to one genomic
        range will be put in a single file with name following
        the pattern chromosome:start-end(strand).
    """
    if not os.path.isdir(directory):
        os.mkdir(directory)
    dump_name = directory + "/" + fname + ".fasta"
    proc = Popen("bedtools getfasta -fi {} -bed stdin -fo stdout -s | fold -w 60 > {}".format(genome, dump_name),
                 stdin=PIPE, shell=True)
    proc.communicate(bed_ranges.encode())
    with open(dump_name, 'r') as infile:
        line = infile.readline()
        while line != "":
            if line.startswith(">"):
                name = line.strip().split(">")[-1]
                with open(directory + "/" + name + '.fasta', 'w') as outfile:
                    outfile.write(line)
                    line = infile.readline()
                    while not line.startswith(">") and line != "":
                        outfile.write(line.upper())
                        line = infile.readline()
            else:
                line = infile.readline()
    return None


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='From given genomes and coordinates get alignments.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument('-leftfile',
                           type=str,
                           nargs='?',
                           dest='leftfilename',
                           help='left genome file')
    argparser.add_argument('-rightfile',
                           type=str,
                           nargs='?',
                           dest='righfilename',
                           help='right genome file')
    argparser.add_argument('-mapfile',
                           type=str,
                           nargs='?',
                           dest='mapfilename',
                           help='map file')
    argparser.add_argument('-outdir',
                           type=str,
                           nargs='?',
                           dest='outdir',
                           help='output directory')

    args = argparser.parse_args()

    left_genome = args.leftfilename
    right_genome = args.righfilename
    synteny_file = args.mapfilename

    with open(synteny_file, 'r') as infile:
        synteny_data = json.load(infile)

    left_rna_bed = get_bed(map(lambda x: dictfilt(x, ['rna_chr', 'rna_start', 'rna_end', 'rna_strand']), synteny_data),
                           'rna_chr', 'rna_start', 'rna_end', 'rna_strand')
    right_synteny_bed = get_bed([item for i in synteny_data for item in i['synteny']])
    get_fasta(right_genome, right_synteny_bed, 'synt', args.outdir)
    get_fasta(left_genome, left_rna_bed, 'rna', args.outdir)
