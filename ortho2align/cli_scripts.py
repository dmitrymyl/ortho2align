import os
import sys
import json
import argparse
import textwrap
import pandas as pd
from pathlib import Path
from collections import defaultdict
from .orthodb import (load_table, filter_table, query_cached_odb_file,
                      split_odb_file, filter_odb_file)


def cache_orthodb_xrefs():
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
                        default='odb10v0_',
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
                      'description']
    genes_dtypes = {'odb_gene_id': 'str',
                    'odb_species_id': 'category',
                    'protein_seq_id': 'str',
                    'synonyms': 'str',
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


def get_orthologs():
    pass


def ortho2align():
    usage = '''\
    ortho2align <command> [options]

    Commands:
        cache_orthodb_xrefs     Cache OrthoDB cross-references into
                                specified directory.

        get_orthodb_map         Retrieve mapping of OrthoDB genes in query
                                species to known orthologs in subject species.

        get_orthologs           Find orthologs of provided query genes in
                                subject species.

    Run ortho2align <command> -h for help on a specific command.
    '''
    commands = {'cache_orthodb_xrefs': cache_orthodb_xrefs,
                'get_orthodb_map': get_orthodb_map,
                'get_orthologs': get_orthologs}
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
