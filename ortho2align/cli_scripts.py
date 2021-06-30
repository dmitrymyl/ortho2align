import argparse
from .pipeline import (bg_from_inter_ranges, bg_from_shuffled_ranges, estimate_background,
                       get_alignments, build_orthologs, get_best_orthologs,
                       annotate_orthologs, run_pipeline)


ortho2align_parser = argparse.ArgumentParser(prog='ortho2align',
                                             description='Ortho2align pipeline',
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                             fromfile_prefix_chars='@')
ortho2align_subparsers = ortho2align_parser.add_subparsers(title='Subcommands',
                                                           metavar='SUBCOMMAND')

bg_from_inter_ranges_parser = ortho2align_subparsers.add_parser('bg_from_inter_ranges',
                                                                help='Generate background set of genomic ranges from intergenic ranges.',
                                                                description='Generate background set of genomic ranges from intergenic ranges.',
                                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
bg_from_inter_ranges_parser.set_defaults(func=bg_from_inter_ranges)
bg_from_inter_ranges_input_group = bg_from_inter_ranges_parser.add_argument_group('Input')
bg_from_inter_ranges_input_group.add_argument('-genes',
                                              type=str,
                                              nargs='?',
                                              required=True,
                                              dest='genes_filename',
                                              help='Gene annotation filename to use for composing intergenic ranges.')
bg_from_inter_ranges_input_group.add_argument('-name_regex',
                                              type=str,
                                              nargs='?',
                                              default=None,
                                              help='Regular expression for extracting gene names from the genes annotation (.gff and .gtf only). '
                                                   'Must contain one catching group.')
bg_from_inter_ranges_processing_group = bg_from_inter_ranges_parser.add_argument_group('Processing')
bg_from_inter_ranges_processing_group.add_argument('-sample_size',
                                                   type=int,
                                                   nargs='?',
                                                   required=True,
                                                   help='Number of background regions to generate.')
bg_from_inter_ranges_processing_group.add_argument('-seed',
                                                   type=int,
                                                   nargs='?',
                                                   default=123,
                                                   help='random seed number for sampling intergenic regions.')
bg_from_inter_ranges_output_group = bg_from_inter_ranges_parser.add_argument_group('Output')
bg_from_inter_ranges_output_group.add_argument('-output',
                                               type=str,
                                               nargs='?',
                                               required=True,
                                               dest='output_filename',
                                               help='Output filename for background regions annotation in bed6 format.')


bg_from_shuffled_ranges_parser = ortho2align_subparsers.add_parser('bg_from_shuffled_ranges',
                                                                   help='Generate background set of genomic ranges from shuffled ranges.',
                                                                   description='Generate background set of genomic ranges from intergenic ranges.',
                                                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
bg_from_shuffled_ranges_parser.set_defaults(func=bg_from_shuffled_ranges)
bg_from_shuffled_ranges_input_group = bg_from_shuffled_ranges_parser.add_argument_group('Input')
bg_from_shuffled_ranges_input_group.add_argument('-genes',
                                                 type=str,
                                                 nargs='?',
                                                 required=True,
                                                 dest='genes_filename',
                                                 help='Gene annotation filename to use for composing shuffled ranges.')
bg_from_shuffled_ranges_input_group.add_argument('-genome',
                                                 type=str,
                                                 nargs='?',
                                                 required=True,
                                                 dest='genome_filename',
                                                 help='Genome file for checking of chromosome sizes.')
bg_from_shuffled_ranges_input_group.add_argument('-name_regex',
                                                 type=str,
                                                 nargs='?',
                                                 default=None,
                                                 help='Regular expression for extracting gene names from the genes annotation (.gff and .gtf only). '
                                                      'Must contain one catching group.')
bg_from_shuffled_ranges_processing_group = bg_from_shuffled_ranges_parser.add_argument_group('Processing')
bg_from_shuffled_ranges_processing_group.add_argument('-sample_size',
                                                      type=int,
                                                      nargs='?',
                                                      required=True,
                                                      help='Number of background regions to generate.')
bg_from_shuffled_ranges_processing_group.add_argument('-seed',
                                                      type=int,
                                                      nargs='?',
                                                      default=123,
                                                      help='random seed number for sampling intergenic regions.')
bg_from_shuffled_ranges_output_group = bg_from_shuffled_ranges_parser.add_argument_group('Output')
bg_from_shuffled_ranges_output_group.add_argument('-output',
                                                  type=str,
                                                  nargs='?',
                                                  required=True,
                                                  dest='output_filename',
                                                  help='Output filename for background regions annotation in bed6 format.')


estimate_background_parser = ortho2align_subparsers.add_parser('estimate_background',
                                                               help='Estimate background alignment scores.',
                                                               description='Estimate background alignment scores.',
                                                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
estimate_background_parser.set_defaults(func=estimate_background)
estimate_background_input_group = estimate_background_parser.add_argument_group('Input')
estimate_background_input_group.add_argument('-query_genes',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             dest='query_genes_filename',
                                             help='Query species gene annotation filename.')
estimate_background_input_group.add_argument('-bg_ranges',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             dest='bg_ranges_filename',
                                             help='Background genomic range set of subject species.')
estimate_background_input_group.add_argument('-query_genome',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             dest='query_genome_filename',
                                             help='Query species genome filename (fasta format).')
estimate_background_input_group.add_argument('-subject_genome',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             dest='subject_genome_filename',
                                             help='Subject species genome filename (fasta format).')
estimate_background_input_group.add_argument('-query_name_regex',
                                             type=str,
                                             nargs='?',
                                             default=None,
                                             help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                                  'Must contain one catching group.')
estimate_background_input_group.add_argument('-bg_name_regex',
                                             type=str,
                                             nargs='?',
                                             default=None,
                                             help='Regular expression for extracting gene names from the background ranges annotation (.gff and .gtf only). '
                                                  'Must contain one catching group.')
estimate_background_parameters_group = estimate_background_parser.add_argument_group('Parameters')
estimate_background_parameters_group.add_argument('-word_size',
                                                  type=int,
                                                  nargs='?',
                                                  default=6,
                                                  help='-word_size argument to use in blastn search.')
estimate_background_parameters_group.add_argument('-observations',
                                                  type=int,
                                                  nargs='?',
                                                  default=1000,
                                                  help='number of scores to retain for each query gene.')
estimate_background_output_group = estimate_background_parser.add_argument_group('Output')
estimate_background_output_group.add_argument('-outdir',
                                              type=str,
                                              nargs='?',
                                              required=True,
                                              help='Output directory name for background files.')
estimate_background_processing_group = estimate_background_parser.add_argument_group('Processing')
estimate_background_processing_group.add_argument('-cores',
                                                  type=int,
                                                  nargs='?',
                                                  default=1,
                                                  help='Number of cores to use for alignment multiprocessing.')
estimate_background_processing_group.add_argument('-seed',
                                                  type=int,
                                                  nargs='?',
                                                  default=123,
                                                  help='random seed for sampling scores.')
estimate_background_processing_group.add_argument('--silent',
                                                  action='store_true',
                                                  help='silent CLI if included.')

get_alignments_parser = ortho2align_subparsers.add_parser('get_alignments',
                                                          help='Compute orthologous alignments of provided query genes and subject species genome.',
                                                          description='Compute orthologous alignments of provided query genes and subject species genome.',
                                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_alignments_parser.set_defaults(func=get_alignments)
get_alignments_input_group = get_alignments_parser.add_argument_group('Input')
get_alignments_input_group.add_argument('-query_genes',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='query genes annotation filename')
get_alignments_input_group.add_argument('-query_genome',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='query species genome filename (fasta)')
get_alignments_input_group.add_argument('-subject_genome',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='subject genome filename (fasta)')
get_alignments_input_group.add_argument('-query_name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
get_alignments_input_group.add_argument('-liftover_chains',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='liftover .chain filename')
get_alignments_params_group = get_alignments_parser.add_argument_group('Parameters')
get_alignments_params_group.add_argument('-min_ratio',
                                         type=float,
                                         nargs='?',
                                         default=0.05,
                                         help='minimal ratio of gene overlapping liftover chain to consider it for liftover')
get_alignments_params_group.add_argument('-word_size',
                                         type=int,
                                         nargs='?',
                                         default=6,
                                         help='-word_size parameter to use in blastn.')
get_alignments_params_group.add_argument('-merge_dist',
                                         type=int,
                                         nargs='?',
                                         default=2000000,
                                         help='how distant two subject syntenies can be to be merged into one syntenic region')
get_alignments_params_group.add_argument('-flank_dist',
                                         type=int,
                                         nargs='?',
                                         default=50000,
                                         help='how many nts to flank syntenic regions in subject species')
get_alignments_output_group = get_alignments_parser.add_argument_group('Output')
get_alignments_output_group.add_argument('-outdir',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='output directory name')
get_alignments_processing_group = get_alignments_parser.add_argument_group('Processing')
get_alignments_processing_group.add_argument('-cores',
                                             type=int,
                                             nargs='?',
                                             default=1,
                                             help='Number of cores to use for alignment multiprocessing')
get_alignments_processing_group.add_argument('--silent',
                                             action='store_true',
                                             help='silent CLI if included.')

build_orthologs_parser = ortho2align_subparsers.add_parser('build_orthologs',
                                                           help='Asses orthologous alignments based on chosen statistical strategy and build orthologs.',
                                                           description='Asses orthologous alignments based on chosen statistical strategy and build orthologs.',
                                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
build_orthologs_parser.set_defaults(func=build_orthologs)
build_orthologs_input_group = build_orthologs_parser.add_argument_group('Input')
build_orthologs_input_group.add_argument('-alignments',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='Alignments dir produced with ortho2align get_alignments.')
build_orthologs_input_group.add_argument('-background',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='background dir produced with ortho2align estimate_background.')
build_orthologs_parameters_group = build_orthologs_parser.add_argument_group('Parameters')
build_orthologs_parameters_group.add_argument('-fitting',
                                              type=str,
                                              nargs='?',
                                              choices=['kde', 'hist'],
                                              default='kde',
                                              help='approach to fit background distribution (kde: KDE, hist: Histogram).')
build_orthologs_parameters_group.add_argument('-threshold',
                                              type=float,
                                              nargs='?',
                                              default=0.05,
                                              help='p-value threshold to filter HSPs by score.')
build_orthologs_parameters_group.add_argument('--fdr',
                                              action='store_true',
                                              help='use FDR correction for HSP scores.')
build_orthologs_output_group = build_orthologs_parser.add_argument_group('Output')
build_orthologs_output_group.add_argument('-outdir',
                                          type=str,
                                          nargs='?',
                                          required=True,
                                          help='output directory.')
build_orthologs_processing_group = build_orthologs_parser.add_argument_group('Processing')
build_orthologs_processing_group.add_argument('-cores',
                                              type=int,
                                              nargs='?',
                                              default=1,
                                              help='Number of cores to use for refinement multiprocessing.')
build_orthologs_processing_group.add_argument('-timeout',
                                              type=int,
                                              nargs='?',
                                              default=None,
                                              help='Time in seconds to terminate a single process of refinement of a single alignment.')
build_orthologs_processing_group.add_argument('--silent',
                                              action='store_true',
                                              help='silent CLI if included.')

get_best_orthologs_parser = ortho2align_subparsers.add_parser('get_best_orthologs',
                                                              help='Select only one ortholog for each query gene based on provided variety of strategies.',
                                                              description='Select only one ortholog for each query gene based on provided variety of strategies.',
                                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_best_orthologs_parser.set_defaults(func=get_best_orthologs)
get_best_orthologs_input_group = get_best_orthologs_parser.add_argument_group('Input')
get_best_orthologs_input_group.add_argument('-query_orthologs',
                                            type=str,
                                            nargs='?',
                                            required=True,
                                            help='query orthologs bed12 file.')
get_best_orthologs_input_group.add_argument('-subject_orthologs',
                                            type=str,
                                            nargs='?',
                                            required=True,
                                            help='subject orthologs bed12 file.')
get_best_orthologs_parameters_group = get_best_orthologs_parser.add_argument_group('Parameters')
get_best_orthologs_parameters_group.add_argument('-value',
                                                 type=str,
                                                 nargs='?',
                                                 choices=['total_length', 'block_count', 'block_length', 'weight'],
                                                 required=True,
                                                 default='block_length',
                                                 help='which value of orthologs to use in case of multiple orthologs.')
get_best_orthologs_parameters_group.add_argument('-function',
                                                 type=str,
                                                 nargs='?',
                                                 choices=['max', 'min'],
                                                 required=True,
                                                 help='orthologs with which value to select in case of multiple orthologs.')
get_best_orthologs_output_group = get_best_orthologs_parser.add_argument_group('Output')
get_best_orthologs_output_group.add_argument('-outfile_query',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             help='output filename for query orthologs.')
get_best_orthologs_output_group.add_argument('-outfile_subject',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             help='output filename for subject orthologs.')
get_best_orthologs_output_group.add_argument('-outfile_map',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             help='output json filename for mapping of query and subject ortholog names.')

annotate_orthologs_parser = ortho2align_subparsers.add_parser('annotate_orthologs',
                                                              help='Annotate found orthologs with provided annotation of subject genome lncRNAs.',
                                                              description='Annotate found orthologs with provided annotation of subject genome lncRNAs.',
                                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
annotate_orthologs_parser.set_defaults(func=annotate_orthologs)
annotate_orthologs_input_group = annotate_orthologs_parser.add_argument_group('Input')
annotate_orthologs_input_group.add_argument('-subject_orthologs',
                                            type=str,
                                            nargs='?',
                                            required=True,
                                            help='subject orthologs filename generated with get_best_orthologs.')
annotate_orthologs_input_group.add_argument('-subject_annotation',
                                            type=str,
                                            nargs='?',
                                            required=True,
                                            help='subject genome lncRNA annotation filename.')
annotate_orthologs_input_group.add_argument('-subject_name_regex',
                                            type=str,
                                            nargs='?',
                                            default=None,
                                            help='Regular expression for extracting gene names from the subject genome lncRNA annotation (.gff and .gtf only). '
                                                 'Must contain one catching group.')
annotate_orthologs_output_group = annotate_orthologs_parser.add_argument_group('Output')
annotate_orthologs_output_group.add_argument('-output',
                                             type=str,
                                             nargs='?',
                                             required=True,
                                             help='output filename.')

run_pipeline_parser = ortho2align_subparsers.add_parser('run_pipeline',
                                                        description='Run the whole pipeline.',
                                                        help='Run the whole pipeline.',
                                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                                        fromfile_prefix_chars='@')
run_pipeline_parser.set_defaults(func=run_pipeline)
run_pipeline_input_group = run_pipeline_parser.add_argument_group('Input')
run_pipeline_input_group.add_argument('-query_genes',
                                      type=str,
                                      nargs='?',
                                      required=True,
                                      help='query genomic ranges.')
run_pipeline_input_group.add_argument('-query_genome',
                                      type=str,
                                      nargs='?',
                                      required=True,
                                      help='Query species genome filename (fasta format).')
run_pipeline_input_group.add_argument('-subject_annotation',
                                      type=str,
                                      nargs='?',
                                      required=True,
                                      help='subject species genome annotation for construction of background genomic ranges.')
run_pipeline_input_group.add_argument('-subject_genome',
                                      type=str,
                                      nargs='?',
                                      required=True,
                                      help='Subject species genome filename (fasta format).')
run_pipeline_input_group.add_argument('-query_name_regex',
                                      type=str,
                                      nargs='?',
                                      default=None,
                                      help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                           'Must contain one catching group.')
run_pipeline_input_group.add_argument('-subject_name_regex',
                                      type=str,
                                      nargs='?',
                                      default=None,
                                      help='Regular expression for extracting gene names from the subject genome annotation (.gff and .gtf only). '
                                           'Must contain one catching group.')
run_pipeline_input_group.add_argument('-liftover_chains',
                                      type=str,
                                      nargs='?',
                                      required=True,
                                      help='liftover .chain filename')
run_pipeline_output_group = run_pipeline_parser.add_argument_group('Output')
run_pipeline_output_group.add_argument('-outdir',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='output directory name.')
run_pipeline_processing_group = run_pipeline_parser.add_argument_group('Processing')
run_pipeline_processing_group.add_argument('-cores',
                                           type=int,
                                           nargs='?',
                                           default=1,
                                           help='Number of cores to use for multiprocessing.')
run_pipeline_processing_group.add_argument('-word_size',
                                           type=int,
                                           nargs='?',
                                           default=6,
                                           help='-word_size argument to use in blastn search.')
run_pipeline_processing_group.add_argument('-seed',
                                           type=int,
                                           nargs='?',
                                           default=0,
                                           help='random seed for sampling procedures.')
run_pipeline_processing_group.add_argument('--silent',
                                           action='store_true',
                                           help='silent CLI if included.')
run_pipeline_processing_group.add_argument('--annotate',
                                           action='store_true',
                                           help='If included, will annotate found orthologs with subject annotation.')
run_pipeline_estimate_background_group = run_pipeline_parser.add_argument_group('Estimating background')
run_pipeline_estimate_background_group.add_argument('-sample_size',
                                                    type=int,
                                                    nargs='?',
                                                    default=200,
                                                    help='Number of background regions to generate.')
run_pipeline_estimate_background_group.add_argument('-observations',
                                                    type=int,
                                                    nargs='?',
                                                    default=1000,
                                                    help='maximum number of background scores to retain for each query gene.')
run_pipeline_get_alignments_group = run_pipeline_parser.add_argument_group('Getting alignments')
run_pipeline_get_alignments_group.add_argument('-min_ratio',
                                               type=float,
                                               nargs='?',
                                               default=0.05,
                                               help='minimal ratio of gene overlapping liftover chain to consider it for liftover')
run_pipeline_get_alignments_group.add_argument('-merge_dist',
                                               type=int,
                                               nargs='?',
                                               default=2000000,
                                               help='how distant two subject syntenies can be to be merged into one syntenic region')
run_pipeline_get_alignments_group.add_argument('-flank_dist',
                                               type=int,
                                               nargs='?',
                                               default=50000,
                                               help='how many nts to flank syntenic regions in subject species')
run_pipeline_build_orthologs_group = run_pipeline_parser.add_argument_group('Building orthologs')
run_pipeline_build_orthologs_group.add_argument('-fitting',
                                                type=str,
                                                nargs='?',
                                                choices=['kde', 'hist'],
                                                default='kde',
                                                help='approach to fit background distribution (kde: KDE, hist: Histogram).')
run_pipeline_build_orthologs_group.add_argument('-threshold',
                                                type=float,
                                                nargs='?',
                                                default=0.05,
                                                help='p-value threshold to filter HSPs by score.')
run_pipeline_build_orthologs_group.add_argument('--fdr',
                                                action='store_true',
                                                help='use FDR correction for HSP scores if included.')
run_pipeline_build_orthologs_group.add_argument('-timeout',
                                                type=int,
                                                nargs='?',
                                                default=None,
                                                help='Time in seconds to terminate a single process of refinement of a single alignment.')
run_pipeline_get_best_orthologs_group = run_pipeline_parser.add_argument_group('Getting best orthologs')
run_pipeline_get_best_orthologs_group.add_argument('-value',
                                                   type=str,
                                                   nargs='?',
                                                   choices=['total_length', 'block_count', 'block_length', 'weight'],
                                                   default='block_length',
                                                   help='which value of orthologs to use in case of multiple orthologs.')
run_pipeline_get_best_orthologs_group.add_argument('-function',
                                                   type=str,
                                                   nargs='?',
                                                   choices=['max', 'min'],
                                                   default='max',
                                                   help='orthologs with which value to select in case of multiple orthologs.')


def ortho2align():
    args = ortho2align_parser.parse_args()
    dict_args = vars(args)
    func = dict_args['func']
    del dict_args['func']
    func(**dict_args)
