import os
import sys
import argparse
import textwrap


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


def cache_orthodb_xrefs():

    from pathlib import Path
    from .orthodb import split_odb_file


    parser = argparse.ArgumentParser(description='Cache OrthoDB gene ID cross-'
                                                 'references into specific '
                                                 'directory.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     prog='ortho2align ' + sys.argv[1])
    parser.add_argument('-orthodb_path',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Path to OrthoDB folder.')
    parser.add_argument('-cache_path',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Path to where cross-references will be cached.')
    parser.add_argument('-orthodb_prefix',
                        type=str,
                        nargs='?',
                        default='odb10v0_',
                        help='Prefix of all OrthoDB files in OrthoDB folder.')
    parser.add_argument('-xref_suffix',
                        type=str,
                        nargs='?',
                        default='gene_xrefs.tab',
                        help='Suffix of OrthoDB cross-references file.')

    args = parser.parse_args(sys.argv[2:])
    orthodb_path = Path(args.orthodb_path)
    cache_path = Path(args.cache_path)
    orthodb_prefix = args.orthodb_prefix
    xref_suffix = args.xref_suffix

    gene_xrefs_colnames = ['odb_gene_id',
                           'external_gene_id',
                           'external_db']

    if not cache_path.exists():
        os.mkdir(cache_path)

    split_odb_file(xref_suffix,
                   orthodb_path,
                   orthodb_prefix,
                   'external_db',
                   gene_xrefs_colnames,
                   cache_path)


def get_orthodb_map():

    import json
    import pandas as pd
    from pathlib import Path
    from collections import defaultdict
    from .orthodb import (load_table, filter_table, query_cached_odb_file,
                          split_odb_file, filter_odb_file)

    external_dbs = ['GOterm',
                    'InterPro',
                    'NCBIproteinGI',
                    'UniProt',
                    'ENSEMBL',
                    'NCBIgid',
                    'NCBIgenename']
    parser = argparse.ArgumentParser(description='Retrieve mapping between '
                                                 'OrthoDB genes in query '
                                                 'speices to know orthologs '
                                                 'in subject species.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-query_genes',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Path to a list of query gene IDs,'
                             'one ID per line.')
    parser.add_argument('-outfile',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Output json filename.')
    parser.add_argument('-query_db',
                        type=str,
                        nargs='?',
                        required=True,
                        choices=external_dbs,
                        help='Database from which query gene IDs '
                             'were taken.')
    parser.add_argument('-subject_db',
                        type=str,
                        nargs='?',
                        required=True,
                        choices=external_dbs,
                        help='Database from which return IDs of '
                             'query genes orthologs in subject '
                             'species.')
    parser.add_argument('-level_taxid',
                        type=str,
                        nargs='?',
                        required=True,
                        help='NCBI taxid of level at which '
                             'orthologous groups will be found.')
    parser.add_argument('-subject_taxids',
                        type=str,
                        nargs='+',
                        required=True,
                        help='List of subject species NCBI '
                             'taxids delimited by a whitespace.')
    parser.add_argument('-orthodb_path',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Path to OrthoDB folder.')
    parser.add_argument('-cache_path',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Path where OrthoDB cross-'
                             'references were cached.')
    parser.add_argument('-orthodb_prefix',
                        type=str,
                        nargs='?',
                        default='odb10v1_',
                        help='Prefix of all OrthoDB files in OrthoDB folder.')
    parser.add_argument('-tmp_path',
                        type=str,
                        nargs='?',
                        default='.tmp/',
                        help='Where to put temporary files.')

    args = parser.parse_args(sys.argv[2:])
    query_genes = Path(args.query_genes)
    output_json_file = Path(args.outfile)
    query_db = args.query_db
    subject_db = args.subject_db
    level_taxid = args.level_taxid
    subject_taxids = args.subject_taxids
    orthodb_path = Path(args.orthodb_path)
    cache_path = Path(args.cache_path)
    orthodb_prefix = args.orthodb_prefix
    tmp_path = Path(args.tmp_path)
    process_id = os.getpid()

    if not tmp_path.exists():
        os.mkdir(tmp_path)
    tmp_proc = tmp_path.joinpath(str(process_id))
    os.mkdir(tmp_proc)

    with open(query_genes, 'r') as infile:
        query_accessions = {line.strip() for line in infile}

    # Load xrefs of query genes from cached file.
    gene_xrefs_colnames = ['odb_gene_id',
                           'external_gene_id',
                           'external_db']
    query_cache_filename = query_cached_odb_file("gene_xrefs.tab",
                                                 orthodb_path,
                                                 orthodb_prefix,
                                                 'external_db',
                                                 query_db,
                                                 gene_xrefs_colnames,
                                                 cache_path)
    query_xref_filename = tmp_proc.joinpath('query_xref.tab')

    query_xref_colnames = ['odb_gene_id',
                           'query_xref_id']
    query_xref_dtypes = {'odb_gene_id': 'str',
                         'query_xref_id': 'str'}

    with open(query_cache_filename, 'r') as infile:
        with open(query_xref_filename, 'w') as outfile:
            filter_table(infile,
                         outfile,
                         'query_xref_id',
                         query_accessions,
                         query_xref_colnames,
                         '\t')
    query_genes_with_odb_id = load_table(query_xref_filename,
                                         query_xref_colnames,
                                         query_xref_dtypes)

    # Leveling OGs.
    OGs_colnames = ['og_id',
                    'level_taxid',
                    'og_name']

    leveled_ogs_filename = query_cached_odb_file('OGs.tab',
                                                 orthodb_path,
                                                 orthodb_prefix,
                                                 'level_taxid',
                                                 level_taxid,
                                                 OGs_colnames,
                                                 tmp_proc)

    leveled_OGs_colnames = ['og_id',
                            'og_name']
    leveled_OGs_dtypes = {'og_id': 'str',
                          'og_name': 'str'}
    leveled_OGs = load_table(leveled_ogs_filename,
                             leveled_OGs_colnames,
                             leveled_OGs_dtypes,
                             usecols=['og_id'])

    # Finding OGs for query genes.
    OG2genes_colnames = ['og_id',
                         'odb_gene_id']
    OG2genes_dtypes = {'og_id': 'category',
                       'odb_gene_id': 'category'}
    genes_in_leveled_OGs_filename = filter_odb_file('OG2genes.tab',
                                                    orthodb_path,
                                                    orthodb_prefix,
                                                    'og_id',
                                                    set(leveled_OGs['og_id']),
                                                    OG2genes_colnames,
                                                    tmp_proc,
                                                    'leveled_OGs')

    OGs_for_query_genes_filename = tmp_proc.joinpath('query_genes_OGs.tab')
    with open(genes_in_leveled_OGs_filename, 'r') as infile:
        with open(OGs_for_query_genes_filename, 'w') as outfile:
            filter_table(infile,
                         outfile,
                         'odb_gene_id',
                         set(query_genes_with_odb_id['odb_gene_id']),
                         OG2genes_colnames,
                         '\t')

    OGs_for_query_genes = load_table(OGs_for_query_genes_filename,
                                     OG2genes_colnames,
                                     OG2genes_dtypes)
    OGs_for_query_genes = pd.merge(OGs_for_query_genes,
                                   query_genes_with_odb_id,
                                   on='odb_gene_id')

    # Getting OrthoDB IDs for subject species.
    species_colnames = ['ncbi_tax_id',
                        'odb_species_id',
                        'name',
                        'genome_assembly_id',
                        'clustered_count',
                        'OG_count',
                        'mapping']
    species_dtypes = {'ncbi_tax_id': 'str',
                      'ob_species_id': 'str',
                      'name': 'str',
                      'genome_assembly_id': 'str',
                      'clustered_count': 'int',
                      'OG_count': 'int',
                      'mapping': 'category'}
    odb_subject_species_ids_filename = filter_odb_file('species.tab',
                                                       orthodb_path,
                                                       orthodb_prefix,
                                                       'ncbi_tax_id',
                                                       subject_taxids,
                                                       species_colnames,
                                                       tmp_proc,
                                                       'species')
    odb_subject_species_ids = load_table(odb_subject_species_ids_filename,
                                         species_colnames,
                                         species_dtypes,
                                         usecols=['ncbi_tax_id',
                                                  'odb_species_id'])

    # Getting subject species genes.
    genes_colnames = ['odb_gene_id',
                      'odb_species_id',
                      'protein_seq_id',
                      'synonyms',
                      'UniProt',
                      'ENSEMBL',
                      'NCBIgid',
                      'description']
    genes_dtypes = {'odb_gene_id': 'str',
                    'odb_species_id': 'category',
                    'protein_seq_id': 'str',
                    'synonyms': 'str',
                    'UniProt': 'str',
                    'ENSEMBL': 'str',
                    'NCBIgid': 'str',
                    'description': 'str'}

    subject_species_genes_filename = filter_odb_file('genes.tab',
                                                     orthodb_path,
                                                     orthodb_prefix,
                                                     'odb_species_id',
                                                     set(odb_subject_species_ids['odb_species_id']),
                                                     genes_colnames,
                                                     tmp_proc,
                                                     'species')
    subject_species_genes = load_table(subject_species_genes_filename,
                                       genes_colnames,
                                       genes_dtypes,
                                       usecols=['odb_gene_id',
                                                'odb_species_id'])

    # Getting subject species genes in leveled OGs.
    subject_genes_in_leveled_OGs_filename = tmp_proc.joinpath('subject_genes_in_leveled_OGs.tab')
    with open(genes_in_leveled_OGs_filename, 'r') as infile:
        with open(subject_genes_in_leveled_OGs_filename, 'w') as outfile:
            filter_table(infile,
                         outfile,
                         'odb_gene_id',
                         set(subject_species_genes['odb_gene_id']),
                         OG2genes_colnames,
                         '\t')
    subject_genes_in_leveled_OGs = load_table(subject_genes_in_leveled_OGs_filename,
                                              OG2genes_colnames,
                                              OG2genes_dtypes)
    subject_genes_in_leveled_OGs = pd.merge(subject_genes_in_leveled_OGs,
                                            subject_species_genes,
                                            on='odb_gene_id')

    # Finding orthologs for query genes.
    orthologs_for_query_genes = pd.merge(OGs_for_query_genes,
                                         subject_genes_in_leveled_OGs,
                                         on='og_id',
                                         suffixes=('', '_ortholog'))

    # Getting odb_gene_id for subject xref db.
    subject_cache_filename = query_cached_odb_file("gene_xrefs.tab",
                                                   orthodb_path,
                                                   orthodb_prefix,
                                                   'external_db',
                                                   subject_db,
                                                   gene_xrefs_colnames,
                                                   cache_path)

    # Filtering subject xref db with found orthologs.
    subject_xref_colnames = ['odb_gene_id',
                             'subject_xref_id']
    subject_xref_dtypes = {'odb_gene_id': 'str',
                           'subject_xref_id': 'str'}
    orthologs_xrefs_filename = tmp_proc.joinpath('orthologs_xrefs.tab')

    with open(subject_cache_filename, 'r') as infile:
        with open(orthologs_xrefs_filename, 'w') as outfile:
            filter_table(infile,
                         outfile,
                         'odb_gene_id',
                         set(orthologs_for_query_genes['odb_gene_id_ortholog']),
                         subject_xref_colnames,
                         '\t')

    orthologs_xrefs = load_table(orthologs_xrefs_filename,
                                 subject_xref_colnames,
                                 subject_xref_dtypes)

    # Finding xref for orthologs of query genes.
    orthologs_with_xrefs = pd.merge(orthologs_for_query_genes,
                                    orthologs_xrefs,
                                    left_on='odb_gene_id_ortholog',
                                    right_on='odb_gene_id',
                                    suffixes=("", "_ortholog"))
    # Saving to file.
    species_dict = {item['odb_species_id']: item['ncbi_tax_id']
                    for item in odb_subject_species_ids.to_dict('record')}
    orthologs_with_xrefs = orthologs_with_xrefs.assign(ncbi_tax_id=lambda x: x['odb_species_id'].map(species_dict))
    report_df = (orthologs_with_xrefs.drop_duplicates(subset=['query_xref_id',
                                                             'subject_xref_id',
                                                             'odb_species_id'])
                                     .drop(['odb_gene_id_ortholog',
                                            'odb_gene_id',
                                            'odb_species_id'],
                                           axis=1))
    map_df = (report_df.drop(['og_id'], axis=1)
                       .groupby(['query_xref_id', 'ncbi_tax_id'])
                       .agg(list))
    json_map = defaultdict(dict)
    for key, value in map_df.to_dict()['subject_xref_id'].items():
        subject_xref_id = key[0]
        ncbi_tax_id = key[1]
        json_map[subject_xref_id][ncbi_tax_id] = value

    with open(output_json_file, 'w') as outfile:
        json.dump(json_map, outfile)

    # Remove temporary files and directories.
    for file in tmp_proc.glob("*"):
        os.remove(file)
    os.rmdir(tmp_proc)


def align_two_ranges_blast(query, subject, **kwargs):
    return query.align_blast(subject, **kwargs)


def estimate_background():
    import random
    import numpy as np
    from .genomicranges import GenomicRangesList
    from .parallel import NonExceptionalProcessPool

    parser = argparse.ArgumentParser(description='Estimate background alignment scores.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-query_genes',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Query species gene annotation filename.')
    parser.add_argument('-subject_genes',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Subject species gene annotation filename.')
    parser.add_argument('-query_genome',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Query species genome filename (fasta format).')
    parser.add_argument('-subject_genome',
                        type=str,
                        nargs='?',
                        required=True,
                        help='Subject species genome filename (fasta format).')
    parser.add_argument('-observations',
                        type=int,
                        nargs='?',
                        required=True,
                        help='Number of alignments to estimate background.')
    parser.add_argument('-output',
                        type=str,
                        nargs='?',
                        required=True,
                        help='File to save background alignment scores. Uses numpy.save function.')
    parser.add_argument('-cores',
                        type=int,
                        nargs='?',
                        default=1,
                        help='Number of cores to use for alignment multiprocessing.')
    parser.add_argument('-seed',
                        type=int,
                        nargs='?',
                        default=123,
                        help='random seed number for sampling query genes and subject intergenic regions.')

    args = parser.parse_args(sys.argv[2:])

    random.seed = args.seed

    query_genes_filename = args.query_genes
    query_genes_filetype = query_genes_filename.split('.')[-1]
    subject_genes_filename = args.subject_genes
    subject_genes_filetype = subject_genes_filename.split('.')[-1]
    query_genome_filename = args.query_genome
    subject_genome_filename = args.subject_genome
    observations = args.observations
    output_filename = args.output
    cores = args.cores

    with open(query_genes_filename, 'r') as infile:
        query_genes = GenomicRangesList.parse_annotation(infile,
                                                         fileformat=query_genes_filetype,
                                                         sequence_file_path=query_genome_filename)
    with open(subject_genes_filename, 'r') as infile:
        subject_genes = GenomicRangesList.parse_annotation(infile,
                                                           fileformat=subject_genes_filetype,
                                                           sequence_file_path=subject_genome_filename)
    subject_genes = subject_genes.inter_ranges()
    query_samples = GenomicRangesList([random.choice(query_genes) for _ in range(observations)],
                                      query_genes.sequence_file_path)
    subject_samples = GenomicRangesList([random.choice(subject_genes) for _ in range(observations)],
                                         subject_genes.sequence_file_path)
    query_samples.get_fasta('query_')
    subject_samples.get_fasta('subject_')
    sample_pairs = zip(query_samples, subject_samples)

    with NonExceptionalProcessPool(cores) as p:
        alignments, exceptions = p.starmap_async(align_two_ranges_blast, sample_pairs)

    if len(exceptions) > 0:
        print('Exceptions occured: ')
        for exception in exceptions:
            print(exception)

    scores = [hsp.score for alignment in alignments for hsp in alignment.HSPs]
    np.save(output_filename, scores, allow_pickle=False)


def align_syntenies(grange):
    grange.relations['neighbours'].get_fasta('neigh_')
    grange.relations['syntenies'].get_fasta('synt_')
    neighbourhood = grange.relations['neighbours'][0]
    alignments = list()
    for synteny in grange.relations['syntenies']:
        alignment = neighbourhood.align_blast(synteny)
        alignment = alignment.to_genomic().cut_coordinates(qleft=grange.start,
                                                           qright=grange.end)
        alignment.srange = grange
        alignments.append(alignment)
    return alignments


def get_alignments():
    import json
    from .parallel import NonExceptionalProcessPool
    from .genomicranges import GenomicRangesList, extract_taxid_mapping

    parser = argparse.ArgumentParser(description='Get alignments',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-query_genes',
                        type=str,
                        nargs='?',
                        required=True,
                        help='query genes annotation filename')
    parser.add_argument('-query_proteins',
                        type=str,
                        nargs='?',
                        required=True,
                        help='query proteins annotation filename')
    parser.add_argument('-query_genome',
                        type=str,
                        nargs='?',
                        required=True,
                        help='query species genome filename (fasta)')
    parser.add_argument('-subject_proteins',
                        type=str,
                        nargs='?',
                        required=True,
                        help='subject proteins annotation filename')
    parser.add_argument('-subject_genome',
                        type=str,
                        nargs='?',
                        required=True,
                        help='subject genome filename (fasta)')
    parser.add_argument('-ortho_map',
                        type=str,
                        nargs='?',
                        required=True,
                        help='orthology map filename')
    parser.add_argument('-subject_taxid',
                        type=str,
                        nargs='?',
                        required=True,
                        help='subject species NCBI taxid')
    parser.add_argument('-output',
                        type=str,
                        nargs='?',
                        required=True,
                        help='output filename')
    parser.add_argument('-cores',
                        type=int,
                        nargs='?',
                        default=1,
                        help='Number of cores to use for alignment multiprocessing.')
    parser.add_argument('-neighbour_dist',
                        type=int,
                        nargs='?',
                        default=0,
                        help='distance to seek protein neighbours of query genes')
    parser.add_argument('-flank_dist',
                        type=int,
                        nargs='?',
                        default=0,
                        help='how many nts to flank syntenic regions in subject species')
    parser.add_argument('-merge_dist',
                        type=int,
                        nargs='?',
                        default=0,
                        help='how distant two subject proteins can be to be merged into one syntenic region')

    args = parser.parse_args(sys.argv[2:])

    query_genes_filename = args.query_genes
    query_genes_filetype = query_genes_filename.split('.')[-1]
    query_genome_filename = args.query_genome
    query_proteins_filename = args.query_proteins
    query_proteins_filetype = query_proteins_filename.split('.')[-1]
    subject_proteins_filename = args.subject_proteins
    subject_proteins_filetype = subject_proteins_filename.split('.')[-1]
    subject_genome_filename = args.subject_genome
    ortho_map_filename = args.ortho_map
    output_filename = args.output
    subject_taxid = args.subject_taxid
    cores = args.cores
    neighbour_dist = args.neighbour_dist
    merge_dist = args.merge_dist
    flank_dist = args.flank_dist

    with open(query_genes_filename, 'r') as infile:
        query_genes = GenomicRangesList.parse_annotation(infile,
                                                         fileformat=query_genes_filetype,
                                                         sequence_file_path=query_genome_filename)
    with open(query_proteins_filename, 'r') as infile:
        query_proteins = GenomicRangesList.parse_annotation(infile,
                                                            fileformat=query_proteins_filetype,
                                                            sequence_file_path=query_genome_filename)
    with open(subject_proteins_filename, 'r') as infile:
        subject_proteins = GenomicRangesList.parse_annotation(infile,
                                                              fileformat=subject_proteins_filetype,
                                                              sequence_file_path=subject_genome_filename)
    with open(ortho_map_filename, 'r') as mapfile:
        ortho_map = json.load(mapfile)

    ortho_map = extract_taxid_mapping(ortho_map, subject_taxid)
    query_proteins.relation_mapping(subject_proteins,
                                    ortho_map,
                                    'orthologs')
    query_genes.get_neighbours(query_proteins,
                               neighbour_dist,
                               'neighbours')
    for query_gene in query_genes:
        query_gene.relations['syntenies'] = GenomicRangesList([],
                                                              genome=subject_genome_filename)

        for neighbour in query_gene.relations['neighbours']:
            query_gene.relations['syntenies'].update(neighbour.relations['orthologs'])

        query_gene.relation['syntenies'] = query_genes.relation['syntenies'] \
                                                      .flank(flank_dist) \
                                                      .merge(merge_dist)
        query_gene.relation['neighbours'] = query_gene.relation['neighbours'] \
                                                      .merge(2 * neighbour_dist
                                                             + query_gene.end
                                                             - query_gene.start)

    with NonExceptionalProcessPool(cores) as p:
        alignments, exceptions = p.map(align_syntenies, query_genes)

    if len(exceptions) > 0:
        print('Exceptions occured:')
        for exception in exceptions:
            print(exception)

    alignments = [alignment.to_dict()
                  for group in alignments
                  for alignment in group]

    with open(output_filename, 'w') as outfile:
        json.dump(outfile, alignments)


def refine_alignments():

    import json
    import numpy as np
    from pathlib import Path
    from .genomicranges import GenomicRangesAlignment
    from .parallel import HistogramFitter, KernelFitter

    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-alignments',
                        type=str,
                        nargs='?',
                        required=True,
                        help='input alignments filename produced with ortho2align get_alignments.')
    parser.add_argument('-background',
                        type=str,
                        nargs='?',
                        required=True,
                        help='background file produced with ortho2align estimate_background.')
    parser.add_argument('-fitting',
                        type=str,
                        nargs='?',
                        choices=['kde', 'hist'],
                        default='kde',
                        help='approach to fit background distribution (kde: KDE, hist: Histogram).')
    parser.add_argument('-threshold',
                        type=float,
                        nargs='?',
                        default=0.05,
                        help='p-value threshold to filter HSPs by score.')
    parser.add_argument('-query_output',
                        type=str,
                        nargs='?',
                        required=True,
                        help='output filename for query transcripts.')
    parser.add_argument('-subject_output',
                        type=str,
                        nargs='?',
                        required=True,
                        help='output filename for subject transcripts.')

    args = parser.parse_args(sys.argv[2:])

    alignments_filename = Path(args.alignments)
    background_filename = Path(args.background)
    fitting_type = args.fitting
    query_output_filename = Path(args.query_output)
    subject_output_filename = Path(args.subject_output)
    pval_treshold = args.threshold

    with open(alignments_filename, 'r') as infile:
        alignment_data = json.load(infile)

    alignment_data = [GenomicRangesAlignment.from_dict(item)
                      for item in alignment_data]

    background_data = np.load(background_filename)

    if fitting_type == 'histogram':
        fitter = HistogramFitter
    elif fitting_type == 'kernel':
        fitter = KernelFitter

    background = fitter(background_data)
    score_threshold = background.isf(pval_treshold)
    query_transcripts = list()
    subject_transcripts = list()

    for alignment in alignment_data:
        alignment.filter_by_socre(score_threshold)
        transcript = alignment.best_transcript()
        record = transcript.to_bed12(mode='str')
        query_transcripts.append(record[0])
        subject_transcripts.append(record[1])

    with open(query_output_filename, 'w') as outfile:
        for record in query_transcripts:
            outfile.write(record + "\n")
    with open(subject_output_filename, 'w') as outfile:
        for record in subject_transcripts:
            outfile.write(record + "\n")


def ortho2align():
    usage = '''\
    ortho2align <command> [options]

    Commands:
        cache_orthodb_xrefs     Cache OrthoDB cross-references into
                                specified directory.

        get_orthodb_map         Retrieve mapping of OrthoDB genes in query
                                species to known orthologs in subject species.

        estimate_background     Estimate background distribution of alignment scores
                                between query genes and subject species intergenic
                                regions.

        get_alignments          Find orthologs of provided query genes in
                                subject species.

        refine_alignments       Refine alignments based on chosen strategy and
                                build final orthologous transcripts.

    Run ortho2align <command> -h for help on a specific command.
    '''
    commands = {'cache_orthodb_xrefs': cache_orthodb_xrefs,
                'get_orthodb_map': get_orthodb_map,
                'estimate_background': estimate_background,
                'get_alignments': get_alignments,
                'refine_alignments': refine_alignments}
    parser = argparse.ArgumentParser(description='ortho2align set of programms.',
                                     usage=textwrap.dedent(usage))
    parser.add_argument('command',
                        help='Subcommand to run.')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands.keys():
        print('Unrecognized command.')
        parser.print_help()
        exit(1)
    commands[args.command]()
