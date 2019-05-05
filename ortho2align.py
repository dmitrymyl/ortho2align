import argparse
import sys
import os
from subprocess import Popen


class CLIVerbose:
    """
    A context manager for verbosity before and after execution
    of some commands.
    """

    def __init__(self, messages, target=sys.stdout):
        """
        Args:
            messages (list): list with 2 items containing start message
                and finish message.
            target (fp): where to print messages (default: sys.stdout).
        """
        self.inmessage, self.outmessage = messages
        self.target = target

    def __enter__(self):
        print(self.inmessage, file=self.target)

    def __exit__(self, exc_type, exc_val, exc_tb):
        print(self.outmessage, file=self.target)


if __name__ == '__main__':
    messages = {'all': ['ortho2align',
                        'The work is done!'],
                'protneigh': ["Estimating protein neighbourhood for given RNAs...",
                              "Protein neighbourhood is estimated."],
                'extract': ["Extracting mapping for given species...",
                            "Extracted mapping for given species from the bulk file."],
                'chromsizes': ["Getting chromsizes from the right genome file...",
                               "Recieved chromsizes."],
                'translation': ["Translating genome annotation to json...",
                                "Annotation is translated."],
                'mapping': ["Mapping RNAs to syntenic regions...",
                            "Syntenic regions are estimated."],
                'sequences': ["Extracting RNA and syntenic sequences...",
                              "Sequences are extracted."],
                'alignments': ["Busy with alignments...",
                               "Completed alignments."]
                }
    parser = argparse.ArgumentParser(description="from given files and parameters produces alignments.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-leftgenome',
                        type=str,
                        nargs='?',
                        help='left genome file')
    parser.add_argument('-rightgenome',
                        type=str,
                        nargs='?',
                        help='right genome file')
    parser.add_argument('-rightannot',
                        type=str,
                        nargs='?',
                        help='right genome annotation (gff|gtf)')
    parser.add_argument('-orthomap',
                        type=str,
                        nargs='?',
                        help='orthology map produced by orthodb_query.py')
    parser.add_argument('-species',
                        type=str,
                        nargs='?',
                        help='species name')
    parser.add_argument('-rnafile',
                        type=str,
                        nargs='?',
                        help='file with RNAs in bed-like format.')
    parser.add_argument('-proteinfile',
                        type='str',
                        nargs='?',
                        help='protein annotation in gff/gtf/bed format.')
    parser.add_argument('-neighbourdistance',
                        type=int,
                        nargs='?',
                        help='neighbourhood radius.')
    parser.add_argument('-rnaheaders',
                        type=str,
                        nargs='+',
                        help='headers of rna file.')
    parser.add_argument('-proteinid',
                        type=str,
                        nargs='?',
                        help='regex pattern to retrieve protein ids.')
    parser.add_argument('-rnaid',
                        type=str,
                        nargs='?',
                        help='column name with unique rna ids in rna file.')
    parser.add_argument('-mergedist',
                        type=str,
                        nargs='?',
                        default=1000,
                        help='max distance between two ranges to be merged. default: 1000')
    parser.add_argument('-flankdist',
                        type=str,
                        nargs='?',
                        default=0,
                        help='distance to extend both starts and ends of merged ranges, default: 0')
    parser.add_argument('-tileSize',
                        type=str,
                        nargs='+',
                        default=11,
                        help='BLAT tile size (default: 11)')
    parser.add_argument('-minIdentity',
                        type=str,
                        nargs='+',
                        default=90,
                        help='minimum sequence identity (in percent) for BLAT (default: 90)')
    parser.add_argument('-directory',
                        type=str,
                        nargs='?',
                        default='fasta',
                        help='directory to store retrieved sequences')
    parser.add_argument('-cores',
                        type=str,
                        nargs='?',
                        default=20,
                        help='number of cores for alignment (default: 20)')
    parser.add_argument('-outfile',
                        type=str,
                        nargs='?',
                        dest='outfilename',
                        help='output json with alignments')
    parser.add_argument('--noneighbour',
                        action='store_true',
                        help='if used, no protein neighbourhood will be estimated')
    parser.add_argument('--noextract',
                        action='store_true',
                        help='if used, no extraction from orthomap will be done')
    parser.add_argument('--nochromsizes',
                        action='store_true',
                        help='if used, no getting chromsizes will be done')
    parser.add_argument('--notranslation',
                        action='store_true',
                        help='if used, no annotation to json will be done')
    parser.add_argument('--nomapping',
                        action='store_true',
                        help='if used, no mapping of syntenic regions will be done')
    parser.add_argument('--nosequences',
                        action='store_true',
                        help='if used, no retrieving sequences will be done')
    parser.add_argument('--noalignment',
                        action='store_true',
                        help='if used, no alignment will be done')

    with CLIVerbose(messages['all']):
        args = parser.parse_args()
        species_name = "_".join(args.species.split())
        source_path = os.path.abspath(os.path.dirname(sys.argv[0]))

        if not args.noneighbour:
            with CLIVerbose(messages['protneigh']):
                proc = Popen("python {}get_neighbourhood.py -rnafile {} -proteinfile {} -distance {} -rnaheaders {} -outfile {} -proteinid {} -rnaid {}".format(args.rnafile,
                                                                                                                                                                args.proteinfile,
                                                                                                                                                                args.neighbourdistance,
                                                                                                                                                                " ".join(args.rnaheaders),
                                                                                                                                                                species_name + "_neighbourhood.json",
                                                                                                                                                                args.proteinid,
                                                                                                                                                                args.rnaid),
                             shell=True)
                proc.communicate()

        if not args.noextract:
            with CLIVerbose(messages['extract']):
                proc = Popen("python {}extract_mapping.py -infile {} -outfile {} -species '{}'".format(source_path,
                                                                                                       args.orthomap,
                                                                                                       species_name + "_orthomap.json",
                                                                                                       args.species),
                             shell=True)
                proc.communicate()

        if not args.nochromsizes:
            with CLIVerbose(messages['chromsizes']):
                proc = Popen("python {}chromsizes_fasta.py -genome {} -outfile {}".format(source_path,
                                                                                          args.rightgenome,
                                                                                          species_name + "_chromsizes.json"),
                             shell=True)
                proc.communicate()

        annot_format = args.rightannot.split(".")[-1]
        if not args.notranslation:
            with CLIVerbose(messages['translation']):
                proc = Popen("python {}annotation2json.py -infile {} -format {} -outfile {}".format(source_path,
                                                                                                    args.rightannot,
                                                                                                    annot_format,
                                                                                                    species_name + "_annot.json"),
                             shell=True)
                proc.communicate()

        if not args.nomapping:
            with CLIVerbose(messages['mapping']):
                proc = Popen("python {}map_synteny.py -leftfile {} -rightfile {} -mapfile {} -outfile {} -mergedist {} -flankdist {} -chromsizes {}".format(source_path,
                                                                                                                                                            species_name + "_neighbourhood.json",
                                                                                                                                                            species_name + "_annot.json",
                                                                                                                                                            species_name + "_orthomap.json",
                                                                                                                                                            species_name + "_synteny.json",
                                                                                                                                                            args.mergedist,
                                                                                                                                                            args.flankdist,
                                                                                                                                                            species_name + "_chromsizes.json"),
                             shell=True)
                proc.communicate()

        if not args.nosequences:
            with CLIVerbose(messages['sequences']):
                proc = Popen("python {}get_fasta.py -leftfile {} -rightfile {} -mapfile {} -outdir {}".format(source_path,
                                                                                                             args.leftgenome,
                                                                                                             args.rightgenome,
                                                                                                             species_name + "_synteny.json",
                                                                                                             args.directory),
                             shell=True)
                proc.communicate()

        if not args.noalignment:
            with CLIVerbose(messages['alignments']):
                proc = Popen("python {}grid_alignment.py -directory {} -mapfile {} -outfile {} -tileSize {} -minIdentity {} -cores {}".format(source_path,
                                                                                                                                              args.directory,
                                                                                                                                              species_name + "_synteny.json",
                                                                                                                                              args.outfilename,
                                                                                                                                              " ".join(args.tileSize),
                                                                                                                                              " ".join(args.minIdentity),
                                                                                                                                              args.cores),
                             shell=True)
                proc.communicate()
