import sys
import argparse
import textwrap
from .pipeline import (cache_orthodb_xrefs, get_orthodb_map,
                       get_orthodb_by_taxid, get_liftover_map,
                       bg_from_inter_ranges, estimate_background,
                       get_alignments, build_orthologs, get_best_orthologs,
                       benchmark_orthologs, annotating_pipeline)


def ortho2align():
    usage = '''\
    ortho2align <command> [options]

    Commands:
        run_pipeline
        cache_orthodb_xrefs     Cache OrthoDB cross-references into
                                specified directory.

        get_orthodb_map         Retrieve mapping of OrthoDB genes in query
                                species to known orthologs in subject species.

        get_orthodb_by_taxid    Retrieve OrthoDB mapping for specific taxid from
                                bulk OrthoDB mapping produced with get_orthodb_map.

        get_liftover_map        Retrieve mapping and annotation of syntenic regions
                                from liftOver chain file.

        bg_from_inter_ranges    Generates set of unique background regions
                                from intergenic ranges of provided set of genes.

        estimate_background     Estimate background distribution of alignment scores
                                between query genes and subject species intergenic
                                regions.

        get_alignments          Compute orthologous alignments of provided query genes
                                and subject species genome.

        build_orthologs         Asses orthologous alignments based on chosen statistical
                                strategy and build orthologs.

        get_best_orthologs      Select only one ortholog for each query gene based
                                on provided variety of strategies.

        benchmark_orthologs     Compare found orthologs against real orthologs and
                                calculate several performance metrics.

    Run ortho2align <command> -h for help on a specific command.
    '''
    commands = {'run_pipeline': annotating_pipeline,
                'cache_orthodb_xrefs': cache_orthodb_xrefs,
                'get_orthodb_map': get_orthodb_map,
                'get_orthodb_by_taxid': get_orthodb_by_taxid,
                'get_liftover_map': get_liftover_map,
                'bg_from_inter_ranges': bg_from_inter_ranges,
                'estimate_background': estimate_background,
                'get_alignments': get_alignments,
                'build_orthologs': build_orthologs,
                'get_best_orthologs': get_best_orthologs,
                'benchmark_orthologs': benchmark_orthologs}

    parser = argparse.ArgumentParser(description='ortho2align set of programms.',
                                     usage=textwrap.dedent(usage))
    parser.add_argument('command',
                        help='Subcommand to run.')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands.keys():
        print('Unrecognized command.')
        parser.print_help()
        exit(1)
    commands[args.command](sys.argv[2:])
