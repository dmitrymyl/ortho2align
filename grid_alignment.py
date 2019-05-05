import json
import os
import re
import argparse
import multiprocessing as mp
from subprocess import Popen
from itertools import product
from math import log


get_coord = lambda x, y: "{}:{}-{}({})".format(x[y[0]], x[y[1]], x[y[2]], x.get(y[3], "."))


def psl_parser(string):
    """
    Parsers psl3 output into dictionary.

    Args:
        string (str): a content of .psl file.

    Returns:
        a list of dictionaries matching psl3 format with one alignment
        per dictionary.
    """
    lines = string.split("\n")
    if not lines[0].startswith('psLayout'):
        raise Exception("The file format is not a psl: incorrect first line.")
    if len(lines) < 5:
        raise Exception("The file format is not a psl: incorrect number of lines.")
    remove_spaces = lambda string: string.replace(" ", "")
    firstline = [remove_spaces(colname) for colname in lines[2].split("\t")]
    secondline = [remove_spaces(colname) for colname in lines[3].split("\t")] + [""] * 3
    header = [start + end
              for start, end in zip(firstline, secondline)]
    k = 5
    matches = list()
    while lines[k] != "":
        match = lines[k].split()
        if len(header) != len(match):
            raise Exception("The file format is not a psl: inconsistent number of columns.")
        match_dict = {key: value for key, value in zip(header, match)}
        matches.append(match_dict)
        k += 1
    return matches


def psl_parser_parser(psl_records, flank=False, score_system='both', **kwargs):
    """
    Transforms psl output: brings coordinates to genomic ones and counts
    blat, blast or user defined alignment scores.

    Args:
        psl_records (list): list of dicts from psl_parser. Sequence names must
            follow the pattern: chr:start-end(strand).
        flank (bool): whether to flank alignment results to the query sizes
            or not (default: False).
        score_system (string): which alignment score system use: 'blast',
            'blat' or 'both' (default: 'both'). In case of 'blast' or 'blat'
            the other score will be set to 0.
        match (int, optional): score if bases match.
        mismatch (int, optional): penalty if bases mismatch.
        gapopen (int, optional): penalty for gap opening.
        gapexted (int, optional): penalty for gap extending per base.

    Returns:
        list of dictionaries, that contains one alignment per dictionary.
        Each alignment contains regions of start and end in both query
        and sequence, alignment blocks sizes, % of identity
        and alignment scores.

    Notes:
        By default, blast and blat score systems differ in their values.
        Blast:
            match: 1
            mismatch: -3
            gapopen: -5
            gapextend: -2
        Blat:
            match: 1
            mismatch: -1
            gapopen: -1
            gapextend: 0
        The user can specify its own scores and penalties so that it will be
        applied to both scoring systems.
    """
    parse_name_pattern = r"(.*):(\d*)-(\d*)\((.)\)"
    name_structure = lambda i: re.search(parse_name_pattern, i).groups()
    alignments = []
    for record in psl_records:
        # initial transformations
        qstructure = name_structure(record["Qname"])
        tstructure = name_structure(record["Tname"])
        true_qstarts = list(map(str, map(lambda i: i + int(qstructure[1]), map(int, record["qStarts"].strip(",").split(",")))))
        true_tstarts = list(map(str, map(lambda i: i + int(tstructure[1]), map(int, record["tStarts"].strip(",").split(",")))))

        # dealing with T coordinates
        if flank:
            if record['strand'] == "+":
                tstart = max(int(record["Tstart"]) - int(record['Qstart']), 0)
                tend = int(record['Tend']) + int(record['Qsize']) - int(record['Qend'])
            else:
                tstart = max(int(record['Tstart']) - int(record['Qsize']) + int(record['Qend']), 0)
                tend = int(record['Tend']) + int(record['Tstart'])

        else:
            tstart = str(int(record["Tstart"]) + int(tstructure[1]))
            tend = str(int(record["Tend"]) + int(tstructure[1]))

        # dealing with Q coordinates
        qstart = str(int(record["Qstart"]) + int(qstructure[1]))
        qend = str(int(record["Qend"]) + int(qstructure[1]))

        # counting score
        blat_counts = [kwargs.get('match', 1),
                       kwargs.get('mismatch', -1),
                       kwargs.get('gapopen', -1),
                       kwargs.get('gapextend', 0)]
        blast_counts = [kwargs.get('match', 1),
                        kwargs.get('mismatch', -3),
                        kwargs.get('gapopen', -5),
                        kwargs.get('gapextend', -2)]
        variables = [(int(record['match']) + int(record['rep.match'])),
                     int(record['mis-match']),
                     (int(record['Qgapcount']) + int(record['Tgapcount'])),
                     (int(record['Qgapbases']) + int(record['Tgapbases']))]
        if score_system == 'blat':
            blast_score = 0
            blat_score = sum([factor * var for factor, var in zip(blat_counts, variables)])
        elif score_system == 'blast':
            blast_score = sum([factor * var for factor, var in zip(blast_counts, variables)])
            blat_score = 0
        elif score_system == 'both':
            blast_score = sum([factor * var for factor, var in zip(blast_counts, variables)])
            blat_score = sum([factor * var for factor, var in zip(blat_counts, variables)])
        else:
            blast_score = 0
            blat_score = 0
        # counting identity as in BLAT
        z = max((int(record['Qend']) - int(record['Qstart'])) - (int(record['Tend']) - int(record['Tstart'])), 0)
        t = int(record['match']) + int(record['mis-match']) + int(record['rep.match'])
        d = 1000 * (1 * int(record['mis-match']) + int(record['Qgapcount']) + round(3 * log(1 + z))) / t
        identity = 100 - 0.1 * d

        # combining dictionary for an output
        keys = ['rna_start',
                'rna_end',
                'synt_start',
                'synt_end',
                'synt_strand',
                'rna_starts',
                'synt_starts',
                'blocksizes',
                'blast_score',
                'blat_score',
                'identity']
        values = [qstart,
                  qend,
                  tstart,
                  tend,
                  record['strand'],
                  true_qstarts,
                  true_tstarts,
                  record['blockSizes'].strip(",").split(","),
                  blast_score,
                  blat_score,
                  identity]
        alignments.append({key: value for key, value in zip(keys, values)})
    return alignments


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='From given genomes and coordinates get alignments.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument('-directory',
                           type=str,
                           nargs='?',
                           dest='directory',
                           help='directory with sequences')
    argparser.add_argument('-mapfile',
                           type=str,
                           nargs='?',
                           dest='mapfilename',
                           help='map file')
    argparser.add_argument('-outfile',
                           type=str,
                           nargs='?',
                           dest='outfilename',
                           help='output file')
    argparser.add_argument('-tileSize',
                           type=int,
                           nargs='+',
                           dest='tileSize',
                           default=11,
                           help='BLAT tile size (default: 11)')
    argparser.add_argument('-minIdentity',
                           type=int,
                           nargs='+',
                           dest='minIdentity',
                           default=90,
                           help='minimum sequence identity (in percent) for BLAT (default: 90)')
    argparser.add_argument('-cores',
                           type=int,
                           nargs='?',
                           dest='cores',
                           default=20,
                           help='number of cores to use for alignment (default: 20)')
    argparser.add_argument('-logging',
                           type=int,
                           nargs='?',
                           dest='logging',
                           default=0,
                           help='if 1, then logs alignments, (default: 0)')
    argparser.add_argument('--no_restrict_introns',
                           action='store_false',
                           help='if included, the mismatch size will be restricted to the size of query sequence.')

    args = argparser.parse_args()

    directory = args.directory
    synteny_file = args.mapfilename
    outfilename = args.outfilename

    with open(synteny_file, 'r') as infile:
        synteny_data = json.load(infile)


    def alignment_func(record, tileSize, minIdentity, directory):
        """
        Main function for multiprocess pool object to perform alignments.

        Args:
            record (dict): a record from synteny map file.
            tileSize (int): tile size to use with blat.
            minIdentity (int): minimal identity to use with blat.
            directory (str): a path to the folder containing sequences.

        Returns:
            modified record dictionary with psl records in 'alignment'.
        """
        rna = get_coord(record, ['rna_chr', 'rna_start', 'rna_end', 'rna_strand'])
        proc_id = os.getpid()
        unique_id = "".join(map(str, [proc_id, minIdentity, tileSize, rna]))
        if args.logging:
            print(proc_id, unique_id)
        for synteny in record['synteny']:
            synt = get_coord(synteny, ['chr', 'start', 'end', 'strand'])
            if args.logging:
                print(unique_id)
            if args.no_restrict_introns:
                maxIntron = int(record['rna_end']) - int(record['rna_start'])
            else:
                maxIntron = 750000
            align = Popen('blat -tileSize={} -minIdentity={} -t=dna -q=rna -maxIntron={} "{}/{}.fasta" "{}/{}.fasta" "out.psl-{}"'.format(tileSize,
                                                                                                                                          minIdentity,
                                                                                                                                          maxIntron,
                                                                                                                                          directory,
                                                                                                                                          synt,
                                                                                                                                          directory,
                                                                                                                                          rna,
                                                                                                                                          unique_id),
                          shell=True)
            align.communicate()
            with open('out.psl-{}'.format(unique_id), 'r') as infile:
                synteny['alignment'] = psl_parser_parser(psl_parser(infile.read()))
            os.remove('out.psl-{}'.format(unique_id))
        return record


    for tileSize, minIdentity in product(args.tileSize, args.minIdentity):
        with mp.Pool(args.cores) as p:
            result = p.starmap(alignment_func, ((i, tileSize, minIdentity, directory) for i in synteny_data))
        with open("-".join([outfilename, 'tile', str(tileSize), 'minIdentity', str(minIdentity)]), 'w') as outfile:
            json.dump(result, outfile)
