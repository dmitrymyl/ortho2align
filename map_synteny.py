import argparse
import json


def merge_ranges(range_list, chromsizes, merge_dist=1000, flank_dist=0, remerge=False):
    """
    Taken list of range dictionaries as range_list,
    merge ranges based on merge distance, then extends merged
    ranges on both sides based on flank distance.
    If remerge is True, launch new merge_ranges function with flank_dist = 0
    and merge_dist = 0.
    Args:
        range_list (list): list of genomic range dictionaries with keys
            'chr', 'start', 'end'.
        chromsizes (dict): dict with chromosome sizes.
        merge_dist (int): distance in nt to merge two ranges
            (default: 1000).
        flank_dist (int): number of nt to flank merged genomic ranges
            (default: 0).
        remerge (bool): whether to merge merged genomic ranges or not
            with merge_dict = 0 and flank_dist = 0 (default: False).
            This might be considered in case merged genomic ranges overlap
            after flanking.
    Returns:
        new range_list.
    """
    if not range_list:
        return []
    range_list.sort(key=lambda i: (i['chr'], int(i['start'])))
    buffer_coord = {'chr': range_list[0]['chr'],
                    'start': range_list[0]['start'],
                    'end': range_list[0]['end']
                    }
    coord_holder = list()
    k = 1
    while k < len(range_list):
        if buffer_coord['chr'] == range_list[k]['chr'] and int(buffer_coord['end']) + merge_dist >= int(range_list[k]['start']):
            buffer_coord['end'] = range_list[k]['end']
        else:
            coord_holder.append(buffer_coord)
            buffer_coord = {'chr': range_list[k]['chr'],
                            'start': range_list[k]['start'],
                            'end': range_list[k]['end']
                            }
        k += 1
    coord_holder.append(buffer_coord)
    for group in coord_holder:
        buff_start = str(int(group['start']) - flank_dist)
        group['start'] = buff_start if int(buff_start) > 0 else 0
        buff_end = str(int(group['end']) + flank_dist)
        group['end'] = buff_end if int(buff_end) < int(chromsizes[group['chr']]) else chromsizes[group['chr']]
    if remerge:
        coord_holder = merge_ranges(coord_holder, chromsizes, merge_dist=0)
    return coord_holder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find orthologs of one genome in the other one and compose syntenic ranges.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-leftfile',
                        type=str,
                        nargs='?',
                        dest='leftfilename',
                        help='json file with rnas and surronding proteins of one genome')
    parser.add_argument('-rightfile',
                        type=str,
                        nargs='?',
                        dest='rightfilename',
                        help='genome annotation of the other genome')
    parser.add_argument('-mapfile',
                        type=str,
                        nargs='?',
                        dest='mapfilename',
                        help='json map file')
    parser.add_argument('-outfile',
                        type=str,
                        nargs='?',
                        dest='outfilename',
                        help='output json file')
    parser.add_argument('-mergedist',
                        type=int,
                        nargs='?',
                        dest='merge_dist',
                        default=1000,
                        help='max distance between two ranges to be merged')
    parser.add_argument('-flankdist',
                        type=int,
                        nargs='?',
                        dest='flank_dist',
                        default=0,
                        help='distance to extend both starts and ends of merged ranges')
    parser.add_argument("-chromsizes",
                        type=str,
                        nargs='?',
                        dest='chromsizesfile',
                        help='json with chromsizes')

    args = parser.parse_args()
    with open(args.leftfilename, 'r') as infile:
        left_genes = json.load(infile)
    with open(args.rightfilename, 'r') as infile:
        right_genes = json.load(infile)
    with open(args.mapfilename, 'r') as infile:
        map_genes = json.load(infile)
    with open(args.chromsizesfile, 'r') as infile:
        chromsizes = json.load(infile)

    for rna in left_genes:
        rna['orthologs'] = list()
        for protein in rna['proteins']:
            rna['orthologs'] += [right_genes.get(gene)
                                 for gene in map_genes.get(protein['prot_code'], [])
                                 if right_genes.get(gene)]
        rna['synteny'] = merge_ranges(rna['orthologs'],
                                      chromsizes,
                                      merge_dist=args.merge_dist,
                                      flank_dist=args.flank_dist,
                                      remerge=True)

    with open(args.outfilename, 'w') as outfile:
        json.dump(left_genes, outfile)
