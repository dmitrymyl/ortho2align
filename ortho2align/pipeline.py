import argparse
import json
import os
import random
from tempfile import TemporaryDirectory

from tqdm import tqdm

from .parsing import parse_annotation
from .fitting import ranking
from .genomicranges import GenomicRangesAlignment
from .parallel import ExceptionLogger

cache_orthodb_xrefs_parser = argparse.ArgumentParser(description='Cache OrthoDB gene ID cross-'
                                                                 'references into specific '
                                                                 'directory.',
                                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                                     prog='ortho2align cache_orthodb_xrefs')
cache_orthodb_xrefs_parser.add_argument('-orthodb_path',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Path to OrthoDB folder.')
cache_orthodb_xrefs_parser.add_argument('-cache_path',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Path to where cross-references will be cached.')
cache_orthodb_xrefs_parser.add_argument('-orthodb_prefix',
                                        type=str,
                                        nargs='?',
                                        default='odb10v0_',
                                        help='Prefix of all OrthoDB files in OrthoDB folder.')
cache_orthodb_xrefs_parser.add_argument('-xref_suffix',
                                        type=str,
                                        nargs='?',
                                        default='gene_xrefs.tab',
                                        help='Suffix of OrthoDB cross-references file.')


def cache_orthodb_xrefs(orthodb_path, cache_path, orthodb_prefix, xref_suffix):

    from pathlib import Path

    from .orthodb import split_odb_file

    orthodb_path = Path(orthodb_path)
    cache_path = Path(cache_path)

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


external_dbs = ['GOterm',
                'InterPro',
                'NCBIproteinGI',
                'UniProt',
                'ENSEMBL',
                'NCBIgid',
                'NCBIgenename']

get_orthodb_map_parser = argparse.ArgumentParser(description='Retrieve mapping between '
                                                             'OrthoDB genes in query '
                                                             'speices to know orthologs '
                                                             'in subject species.',
                                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_orthodb_map_parser.add_argument('-query_genes',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='Path to a list of query gene IDs,'
                                         'one ID per line.')
get_orthodb_map_parser.add_argument('-outfile',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='Output json filename.')
get_orthodb_map_parser.add_argument('-query_db',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    choices=external_dbs,
                                    help='Database from which query gene IDs '
                                         'were taken.')
get_orthodb_map_parser.add_argument('-subject_db',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    choices=external_dbs,
                                    help='Database from which return IDs of '
                                         'query genes orthologs in subject '
                                         'species.')
get_orthodb_map_parser.add_argument('-level_taxid',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='NCBI taxid of level at which '
                                         'orthologous groups will be found.')
get_orthodb_map_parser.add_argument('-subject_taxids',
                                    type=str,
                                    nargs='+',
                                    required=True,
                                    help='List of subject species NCBI '
                                         'taxids delimited by a whitespace.')
get_orthodb_map_parser.add_argument('-orthodb_path',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='Path to OrthoDB folder.')
get_orthodb_map_parser.add_argument('-cache_path',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='Path where OrthoDB cross-'
                                         'references were cached.')
get_orthodb_map_parser.add_argument('-orthodb_prefix',
                                    type=str,
                                    nargs='?',
                                    default='odb10v1_',
                                    help='Prefix of all OrthoDB files in OrthoDB folder.')
get_orthodb_map_parser.add_argument('-tmp_path',
                                    type=str,
                                    nargs='?',
                                    default='.tmp/',
                                    help='Where to put temporary files.')
get_orthodb_map_parser.add_argument('--silent',
                                    action='store_true',
                                    help='silent CLI if included.')


def get_orthodb_map(query_genes, output_json_file, query_db, subject_db,
                    level_taxid, subject_taxids, orthodb_path, cache_path,
                    orthodb_prefix, silent, tmp_path):

    import json
    from collections import defaultdict
    from pathlib import Path

    import pandas as pd

    from .orthodb import (filter_odb_file, filter_table, load_table,
                          query_cached_odb_file)

    query_genes = Path(query_genes)
    output_json_file = Path(output_json_file)
    orthodb_path = Path(orthodb_path)
    cache_path = Path(cache_path)
    tmp_path = Path(tmp_path)
    process_id = os.getpid()

    if not tmp_path.exists():
        os.mkdir(tmp_path)
    tmp_proc = tmp_path.joinpath(str(process_id))
    os.mkdir(tmp_proc)

    with open(query_genes, 'r') as infile:
        query_accessions = {line.strip() for line in infile}

    cmd_hints = ['loading xrefs of query genes from cached file...',
                 'levelling OGs...',
                 'finding OGs for query genes...',
                 'getting OrthoDB IDs for subject species...',
                 'getting subject species genes...',
                 'getting subject species genes in leveled OGs...',
                 'finding orthologs for query genes...',
                 'getting odb_gene_id for subject xref db...',
                 'filtering subject xref db with found orthologs...',
                 'finding xref for orthologs of query genes...',
                 'saving to file...',
                 'removing temporary files and directories...',
                 'finished']
    cmd_point = 0
    with tqdm(total=(len(cmd_hints) - 1),
              bar_format='{n_fmt}/{total_fmt} {elapsed}<{remaining} {postfix}',
              postfix=cmd_hints[cmd_point],
              disable=silent) as pbar:

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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
        orthologs_for_query_genes = pd.merge(OGs_for_query_genes,
                                             subject_genes_in_leveled_OGs,
                                             on='og_id',
                                             suffixes=('', '_ortholog'))

        # Getting odb_gene_id for subject xref db.
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
        subject_cache_filename = query_cached_odb_file("gene_xrefs.tab",
                                                       orthodb_path,
                                                       orthodb_prefix,
                                                       'external_db',
                                                       subject_db,
                                                       gene_xrefs_colnames,
                                                       cache_path)

        # Filtering subject xref db with found orthologs.
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
        orthologs_with_xrefs = pd.merge(orthologs_for_query_genes,
                                        orthologs_xrefs,
                                        left_on='odb_gene_id_ortholog',
                                        right_on='odb_gene_id',
                                        suffixes=("", "_ortholog"))
        # Saving to file.
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
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
                           .agg(list)
                           .applymap(lambda x: x if isinstance(x, list) else []))
        json_map = defaultdict(dict)
        for key, value in map_df.to_dict()['subject_xref_id'].items():
            subject_xref_id = key[0]
            ncbi_tax_id = key[1]
            json_map[subject_xref_id][ncbi_tax_id] = value

        with open(output_json_file, 'w') as outfile:
            json.dump(json_map, outfile)

        # Remove temporary files and directories.
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()
        for file in tmp_proc.glob("*"):
            os.remove(file)
        os.rmdir(tmp_proc)
        # Finished.
        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()


get_orthodb_by_taxid_parser = argparse.ArgumentParser(description='Retrieve OrthoDB mapping for specific taxid from bulk OrthoDB mapping produced with get_orthodb_map.',
                                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_orthodb_by_taxid_parser.add_argument('-orthodb_map',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='OrthoDB mapping produced with get_orthodb_map')
get_orthodb_by_taxid_parser.add_argument('-taxid',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='NCBI taxid to extract mapping for')
get_orthodb_by_taxid_parser.add_argument('-output',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='json output filename')


def get_orthodb_by_taxid(orthodb_map_filename, taxid, output_filename):
    import json
    from pathlib import Path

    from .genomicranges import extract_taxid_mapping

    orthodb_map_filename = Path(orthodb_map_filename)
    output_filename = Path(output_filename)

    with open(orthodb_map_filename, 'r') as infile:
        orthodb_map = json.load(infile)

    taxid_map = extract_taxid_mapping(orthodb_map, taxid)

    with open(output_filename, 'w') as outfile:
        json.dump(taxid_map, outfile)


get_liftover_map_parser = argparse.ArgumentParser(description='Retrieve mapping and annotation of syntenic regions from liftOver chain file.',
                                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_liftover_map_parser.add_argument('-chain_file',
                                     type=str,
                                     nargs='?',
                                     required=True,
                                     help='Input liftOver chain file.')
get_liftover_map_parser.add_argument('-liftover_map',
                                     type=str,
                                     nargs='?',
                                     required=True,
                                     help='Output json map file.')
get_liftover_map_parser.add_argument('-query_anchors',
                                     type=str,
                                     nargs='?',
                                     required=True,
                                     help='output query anchors bed file.')
get_liftover_map_parser.add_argument('-subject_anchors',
                                     type=str,
                                     nargs='?',
                                     required=True,
                                     help='output subject anchors bed file.')


def get_liftover_map(chain_file, liftover_map, query_anchors, subject_anchors):
    import json
    from collections import defaultdict, namedtuple
    from pathlib import Path

    chain_filename = Path(chain_file)
    liftover_map_filename = Path(liftover_map)
    query_anchors_filename = Path(query_anchors)
    subject_anchors_filename = Path(subject_anchors)

    Chain = namedtuple('Chain',
                       'chain score qchrom qsize qstrand qstart qend schrom ssize sstrand sstart send name')
    QueryBed = namedtuple('QueryBed',
                          'qchrom qstart qend name score qstrand')
    SubjectBed = namedtuple('SubjectBed',
                            'schrom sstart send name score sstrand')
    chain_map = defaultdict(list)

    with open(chain_filename, 'r') as chainfile, \
         open(query_anchors_filename, 'w') as qoutfile, \
         open(subject_anchors_filename, 'w') as soutfile:
        for line in chainfile:
            if line.startswith('chain'):
                chain = Chain(*line.strip().split())
                query_bed = QueryBed(*(chain._asdict()[k] for k in QueryBed._fields))
                subject_bed = SubjectBed(*(chain._asdict()[k] for k in SubjectBed._fields))
                chain_map[chain.name].append(chain.name)
                qoutfile.write('\t'.join(query_bed) + '\n')
                soutfile.write('\t'.join(subject_bed) + '\n')

    with open(liftover_map_filename, 'w') as outfile:
        json.dump(chain_map, outfile)


bg_from_inter_ranges_parser = argparse.ArgumentParser(description='Generate background set of genomic ranges from intergenic ranges.',
                                                      prog='ortho2align bg_from_inter_ranges',
                                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

bg_from_inter_ranges_parser.add_argument('-genes',
                                         type=str,
                                         nargs='?',
                                         required=True,
                                         help='Gene annotation filename to use for composing intergenic ranges.')
bg_from_inter_ranges_parser.add_argument('-name_regex',
                                         type=str,
                                         nargs='?',
                                         default=None,
                                         help='Regular expression for extracting gene names from the genes annotation (.gff and .gtf only). '
                                              'Must contain one catching group.')
bg_from_inter_ranges_parser.add_argument('-observations',
                                         type=int,
                                         nargs='?',
                                         required=True,
                                         help='Number of background regions to generate.')
bg_from_inter_ranges_parser.add_argument('-output',
                                         type=str,
                                         nargs='?',
                                         help='Output filename for background regions annotation in bed6 format.')
bg_from_inter_ranges_parser.add_argument('-seed',
                                         type=int,
                                         nargs='?',
                                         default=123,
                                         help='random seed number for sampling intergenic regions.')


def bg_from_inter_ranges(genes_filename, name_regex, observations, output_filename, seed):
    import random

    from .genomicranges import BaseGenomicRangesList
    from .parsing import parse_annotation

    random.seed(seed)

    with open(genes_filename, 'r') as infile:
        genes = parse_annotation(infile, name_regex=name_regex)

    inter_genes = genes.inter_ranges()
    if len(inter_genes) < observations:
        raise ValueError(f'The number of observations ({observations}) '
                         'must be less than the number of intergenic '
                         f'regions ({len(inter_genes)}) derived from genes')
    samples = BaseGenomicRangesList(random.sample(inter_genes, k=observations))

    with open(output_filename, 'w') as outfile:
        samples.to_bed6(outfile)


def _align_two_ranges_blast(query, subject, **kwargs):
    try:
        return query.align_blast(subject, **kwargs)
    except Exception as e:
        raise ExceptionLogger(e, {'query': query,
                                  'subject': subject})


def _estimate_bg_for_single_query(query, bg_ranges, word_size, output_name, sample_size=1000, seed=0):
    all_scores = list()
    for bg_range in bg_ranges:
        alignment = _align_two_ranges_blast(query, bg_range, word_size=word_size)
        scores = [hsp.score for hsp in alignment.HSPs]
        all_scores += scores
    score_size = len(all_scores)
    if score_size > sample_size:
        random.seed(seed)
        all_scores = random.sample(all_scores, sample_size)
    with open(output_name, 'w') as outfile:
        json.dump(all_scores, outfile)
    return output_name, score_size


estimate_background_parser = argparse.ArgumentParser(description='Estimate background alignment scores.',
                                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

estimate_background_parser.add_argument('-query_genes',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Query species gene annotation filename.')
estimate_background_parser.add_argument('-bg_ranges',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Background genomic range set of subject species.')
estimate_background_parser.add_argument('-query_genome',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Query species genome filename (fasta format).')
estimate_background_parser.add_argument('-subject_genome',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Subject species genome filename (fasta format).')
estimate_background_parser.add_argument('-outdir',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='Output directory name for background files.')
estimate_background_parser.add_argument('-cores',
                                        type=int,
                                        nargs='?',
                                        default=1,
                                        help='Number of cores to use for alignment multiprocessing.')
estimate_background_parser.add_argument('-name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
estimate_background_parser.add_argument('-word_size',
                                        type=int,
                                        nargs='?',
                                        default=6,
                                        help='-word_size argument to use in blastn search.')
estimate_background_parser.add_argument('-sample_size',
                                        type=int,
                                        nargs='?',
                                        default=1000,
                                        help='number of scores to retain for each query gene.')
estimate_background_parser.add_argument('-seed',
                                        type=int,
                                        nargs='?',
                                        default=123,
                                        help='random seed for sampling scores.')
estimate_background_parser.add_argument('--silent',
                                        action='store_true',
                                        help='silent CLI if included.')


def estimate_background(query_genes, bg_ranges, query_genome,
                        subject_genome, outdir, cores, name_regex,
                        word_size, sample_size, seed, silent=False):
    from tempfile import TemporaryDirectory

    from .genomicranges import FastaSeqFile, SequencePath
    from .parallel import NonExceptionalProcessPool
    from .parsing import parse_annotation

    query_genes_filename = query_genes
    bg_ranges_filename = bg_ranges
    query_genome_filename = query_genome
    subject_genome_filename = subject_genome

    cmd_hints = ['reading annotations...',
                 'getting sample sequences...',
                 'aligning samples...',
                 'finished.']
    with tqdm(cmd_hints,
              bar_format="{r_bar}",
              postfix=cmd_hints[0],
              initial=1,
              disable=silent) as pbar:

        with open(query_genes_filename, 'r') as infile:
            query_genes = parse_annotation(infile,
                                           sequence_file_path=query_genome_filename,
                                           name_regex=name_regex)
        with open(bg_ranges_filename, 'r') as infile:
            bg_ranges = parse_annotation(infile,
                                         sequence_file_path=subject_genome_filename,
                                         name_regex=name_regex)

        pbar.postfix = pbar.iterable[pbar.n]
        pbar.update()

        query_chromsizes = FastaSeqFile(query_genome_filename).chromsizes
        subject_chromsizes = FastaSeqFile(subject_genome_filename).chromsizes

        with TemporaryDirectory() as tempdirname:
            tempdir = SequencePath(tempdirname)

            query_genes.get_fasta(outfileprefix='query_',
                                  outdir=tempdir,
                                  chromsizes=query_chromsizes)
            bg_ranges.get_fasta(outfileprefix='bg_',
                                outdir=tempdir,
                                chromsizes=subject_chromsizes)

            pbar.postfix = pbar.iterable[pbar.n]
            pbar.update()

            if not os.path.exists(outdir):
                os.mkdir(outdir)
            data = ((query,
                     bg_ranges,
                     word_size,
                     os.path.join(outdir, f"{query.name}.json"),
                     sample_size,
                     seed)
                    for query in query_genes)
            try:
                with NonExceptionalProcessPool(cores, verbose=(not silent)) as p:
                    results, exceptions = p.starmap_async(_estimate_bg_for_single_query, data)
                if len(exceptions) > 0:
                    tqdm.write(f'{len(exceptions)} exceptions occured.')
                    for exception in exceptions:
                        tqdm.write(str(exception))
            except Exception as e:
                raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')

        pbar.postfix = pbar.iterable[pbar.n]
        pbar.update()


def _anchor(query_genes, query_anchors, query_genome_filename,
            subject_anchors, subject_genome_filename, subject_chromsizes,
            ortho_map, neighbour_dist, merge_dist, flank_dist):
    from .genomicranges import GenomicRangesList

    query_anchors.relation_mapping(subject_anchors,
                                   ortho_map,
                                   'orthologs')
    query_genes.get_neighbours(query_anchors,
                               neighbour_dist,
                               'neighbours')
    query_unalignable = list()
    query_prepared = list()

    for query_gene in tqdm(query_genes):
        if len(query_gene.relations['neighbours']) == 0:
            query_unalignable.append(query_gene)
            continue
        syntenies = [grange
                     for neighbour in query_gene.relations['neighbours']
                     for grange in neighbour.relations['orthologs']]
        synteny_list = GenomicRangesList(syntenies,
                                         sequence_file_path=subject_genome_filename)
        query_gene.relations['syntenies'] = synteny_list.close_merge(merge_dist).flank(flank_dist,
                                                                                       chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = query_gene.relations['neighbours'].close_merge(float('inf'))
        neighbourhood = query_gene.relations['neighbours'][0]
        if query_gene.end > neighbourhood.end or query_gene.start < neighbourhood.start:
            query_gene.relations['syntenies'] = query_gene.relations['syntenies'].flank(neighbour_dist + query_gene.end - query_gene.start,
                                                                                        chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = GenomicRangesList([query_gene.flank(flank_dist)],
                                                               sequence_file_path=query_genome_filename)

        query_prepared.append(query_gene)

    query_prepared_genes = GenomicRangesList(query_prepared,
                                             sequence_file_path=query_genome_filename)
    query_unalignable_genes = GenomicRangesList(query_unalignable,
                                                sequence_file_path=query_genome_filename)
    return query_prepared_genes, query_unalignable_genes


def _lift(query_genes, query_genome_filename, subject_genome_filename,
          subject_chromsizes, liftover_chains,
          merge_dist, flank_dist, min_ratio=0.05, allow_duplications=True):
    from subprocess import run

    from .genomicranges import GenomicRangesList

    query_unalignable = list()
    query_prepared = list()

    with TemporaryDirectory() as tempdirname:
        temp_query_filename = os.path.join(tempdirname, 'granges.bed')
        temp_result_filename = os.path.join(tempdirname, 'result.bed')
        temp_unmapped_filename = os.path.join(tempdirname, 'unmapped.bed')
        with open(temp_query_filename, 'w') as outfile:
            query_genes.to_bed6(outfile)
        run(f"liftOver -minMatch={min_ratio} -multiple -noSerial {temp_query_filename} {liftover_chains} {temp_result_filename} {temp_unmapped_filename}", shell=True)
        with open(temp_result_filename, 'r') as infile:
            lifted = parse_annotation(infile)

    for query_gene in tqdm(query_genes):
        if lifted.name_mapping.get(query_gene.name) is None:
            query_unalignable.append(query_gene)
            continue
        syntenies = [grange
                     for grange in lifted.name_mapping.get(query_gene.name)]
        synteny_list = GenomicRangesList(syntenies,
                                         sequence_file_path=subject_genome_filename)
        query_gene.relations['syntenies'] = synteny_list.close_merge(merge_dist).flank(flank_dist,
                                                                                       chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = GenomicRangesList([query_gene.flank(flank_dist)],
                                                               sequence_file_path=query_genome_filename)
        query_prepared.append(query_gene)

    query_prepared_genes = GenomicRangesList(query_prepared,
                                             sequence_file_path=query_genome_filename)
    query_unalignable_genes = GenomicRangesList(query_unalignable,
                                                sequence_file_path=query_genome_filename)
    return query_prepared_genes, query_unalignable_genes


def _align_syntenies(grange, word_size, outdir):
    try:
        neighbourhood = grange.relations['neighbours'][0]
        alignments = list()
        for synteny in grange.relations['syntenies']:
            alignment = neighbourhood.align_blast(synteny, word_size=word_size)
            alignment.to_genomic()
            alignment = alignment.cut_coordinates(qleft=grange.start,
                                                  qright=grange.end)
            alignment.qrange = grange
            alignments.append(alignment.to_dict())

        with open(os.path.join(outdir, f"{grange.name}.json"), 'w') as outfile:
            json.dump(alignments, outfile)

        return grange.name, len(alignments)

    except Exception as e:
        raise ExceptionLogger(e, grange)


get_alignments_parser = argparse.ArgumentParser(description='Compute orthologous alignments of provided query genes and subject species genome.',
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_alignments_parser.add_argument('-mode',
                                   type=str,
                                   nargs='?',
                                   choices=['anchor', 'lift'],
                                   default='lift',
                                   help='type of orthology map and neighbourhood size.')
get_alignments_parser.add_argument('-query_genes',
                                   type=str,
                                   nargs='?',
                                   required=True,
                                   help='query genes annotation filename')
get_alignments_parser.add_argument('-query_name_regex',
                                   type=str,
                                   nargs='?',
                                   default=None,
                                   help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                        'Must contain one catching group.')
get_alignments_parser.add_argument('-query_anchors',
                                   type=str,
                                   nargs='?',
                                   help='query anchors annotation filename')
get_alignments_parser.add_argument('-query_anchors_name_regex',
                                   type=str,
                                   nargs='?',
                                   default=None,
                                   help='Regular expression for extracting gene names from the query anchors annotation (.gff and .gtf only). '
                                        'Must contain one catching group.')
get_alignments_parser.add_argument('-query_genome',
                                   type=str,
                                   nargs='?',
                                   required=True,
                                   help='query species genome filename (fasta)')
get_alignments_parser.add_argument('-subject_anchors',
                                   type=str,
                                   nargs='?',
                                   help='subject anchors annotation filename')
get_alignments_parser.add_argument('-subject_anchors_name_regex',
                                   type=str,
                                   nargs='?',
                                   default=None,
                                   help='Regular expression for extracting gene names from the subject anchors annotation (.gff and .gtf only). '
                                        'Must contain one catching group.')
get_alignments_parser.add_argument('-subject_genome',
                                   type=str,
                                   nargs='?',
                                   required=True,
                                   help='subject genome filename (fasta)')
get_alignments_parser.add_argument('-ortho_map',
                                   type=str,
                                   nargs='?',
                                   help='orthology map filename')
get_alignments_parser.add_argument('-liftover_chains',
                                   type=str,
                                   nargs='?',
                                   help='liftover .chain filename')
get_alignments_parser.add_argument('-outdir',
                                   type=str,
                                   nargs='?',
                                   required=True,
                                   help='output directory name')
get_alignments_parser.add_argument('-cores',
                                   type=int,
                                   nargs='?',
                                   default=1,
                                   help='Number of cores to use for alignment multiprocessing')
get_alignments_parser.add_argument('-min_ratio',
                                   type=float,
                                   nargs='?',
                                   default=0.05,
                                   help='minimal ratio of gene overlapping liftover chain to consider it for liftover')
get_alignments_parser.add_argument('-neighbour_dist',
                                   type=int,
                                   nargs='?',
                                   default=0,
                                   help='distance to seek anchor neighbours of query genes')
get_alignments_parser.add_argument('-flank_dist',
                                   type=int,
                                   nargs='?',
                                   default=0,
                                   help='how many nts to flank syntenic regions in subject species')
get_alignments_parser.add_argument('-merge_dist',
                                   type=int,
                                   nargs='?',
                                   default=0,
                                   help='how distant two subject anchors can be to be merged into one syntenic region')
get_alignments_parser.add_argument('-word_size',
                                   type=int,
                                   nargs='?',
                                   default=6,
                                   help='-word_size parameter to use in blastn.')
get_alignments_parser.add_argument('--silent',
                                   action='store_true',
                                   help='silent CLI if included.')


def get_alignments(mode,
                   query_genes,
                   query_genome,
                   query_anchors,
                   subject_anchors,
                   subject_genome,
                   ortho_map,
                   liftover_chains,
                   outdir,
                   query_name_regex=None,
                   query_anchors_name_regex=None,
                   subject_anchors_name_regex=None,
                   cores=1,
                   min_ratio=0.05,
                   neighbour_dist=0,
                   merge_dist=0,
                   flank_dist=0,
                   word_size=6,
                   silent=False):
    import json
    from functools import partial
    from tempfile import TemporaryDirectory

    from .genomicranges import (BaseGenomicRangesList, FastaSeqFile,
                                SequencePath)
    from .parallel import NonExceptionalProcessPool
    from .parsing import parse_annotation

    program_mode = mode
    query_genes_filename = query_genes
    query_genome_filename = query_genome
    query_anchors_filename = query_anchors
    subject_anchors_filename = subject_anchors
    subject_genome_filename = subject_genome
    ortho_map_filename = ortho_map
    liftover_chains_filename = liftover_chains

    # if program_mode == 'lift':
    #     if liftover_chains_filename is None:
    #         parser.error(f'You must supply -liftover_chains with -mode {program_mode}.')

    # if program_mode == 'anchor':
    #     if query_anchors_filename is None:
    #         parser.error(f'You must supply -query_anchors with -mode {program_mode}.')
    #     if subject_anchors_filename is None:
    #         parser.error(f'You must supply -subject_anchors with -mode {program_mode}.')
    #     if ortho_map_filename is None:
    #         parser.error(f'You must supply -ortho_map with -mode {program_mode}.')

    cmd_hints = ['reading annotations...',
                 'mapping orthology and neighbour relations, composing syntenies...',
                 'getting sequences...',
                 'getting alignments...',
                 'saving alignments...',
                 'finished.']
    cmd_point = 0
    with tqdm(total=(len(cmd_hints) - 1),
              bar_format='{n_fmt}/{total_fmt} {elapsed}<{remaining} {postfix}',
              postfix=cmd_hints[cmd_point],
              disable=silent) as pbar:

        with open(query_genes_filename, 'r') as infile:
            query_genes = parse_annotation(infile,
                                           sequence_file_path=query_genome_filename,
                                           name_regex=query_name_regex)

        subject_genome = FastaSeqFile(subject_genome_filename)
        subject_chromsizes = subject_genome.chromsizes

        if program_mode == "anchor":
            with open(query_anchors_filename, 'r') as infile:
                query_anchors = parse_annotation(infile,
                                                 sequence_file_path=query_genome_filename,
                                                 name_regex=query_anchors_name_regex)

            with open(subject_anchors_filename, 'r') as infile:
                subject_anchors = parse_annotation(infile,
                                                   sequence_file_path=subject_genome_filename,
                                                   name_regex=subject_anchors_name_regex)

            with open(ortho_map_filename, 'r') as mapfile:
                ortho_map = json.load(mapfile)

            cmd_point += 1
            pbar.postfix = cmd_hints[cmd_point]
            pbar.update()

            process_args = [query_genes, query_anchors, query_genome_filename,
                            subject_anchors, subject_genome_filename, subject_chromsizes,
                            ortho_map, neighbour_dist, merge_dist, flank_dist]

            query_prepared_genes, query_unalignable_genes = _anchor(*process_args)

        elif program_mode == "lift":
            cmd_point += 1
            pbar.postfix = cmd_hints[cmd_point]
            pbar.update()

            process_args = [query_genes, query_genome_filename, subject_genome_filename,
                            subject_chromsizes, liftover_chains_filename,
                            merge_dist, flank_dist, min_ratio]
            query_prepared_genes, query_unalignable_genes = _lift(*process_args)
        else:
            raise ValueError("program_mode is not one of ['anchor', 'lift']")

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()

        query_chromsizes = query_prepared_genes.sequence_file.chromsizes

        with TemporaryDirectory() as tempdirname:
            tempdir = SequencePath(tempdirname)
            for query_gene in tqdm(query_prepared_genes):
                query_gene.relations['neighbours'].get_fasta(outfileprefix='neigh_',
                                                             outdir=tempdir,
                                                             chromsizes=query_chromsizes)
                query_gene.relations['syntenies'].get_fasta(outfileprefix='synt_',
                                                            outdir=tempdir,
                                                            chromsizes=subject_chromsizes)

            cmd_point += 1
            pbar.postfix = cmd_hints[cmd_point]
            pbar.update()

            if not os.path.exists(outdir):
                os.mkdir(outdir)

            align_syntenies_word_size = partial(_align_syntenies, word_size=word_size, outdir=outdir)

            with NonExceptionalProcessPool(cores, verbose=not silent) as p:
                try:
                    results, exceptions = p.map_async(align_syntenies_word_size, query_prepared_genes)
                except Exception as e:
                    raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')

        query_exception_ranges = list()
        if len(exceptions) > 0:
            tqdm.write(f'{len(exceptions)} exceptions occured:')
            for exception in exceptions:
                tqdm.write(str(exception))
                query_exception_ranges.append(exception.variable)

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()

        query_exception_list = BaseGenomicRangesList(query_exception_ranges)

        with open(os.path.join(outdir, "unalignable.bed"), 'w') as outfile:
            query_unalignable_genes.to_bed6(outfile)

        with open(os.path.join(outdir, "exceptions.bed"), 'w') as outfile:
            query_exception_list.to_bed6(outfile)

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()


def _build(query_gene_alignments_filename, query_bg_filename, fitter, fdr, pval_threshold):
    with open(query_gene_alignments_filename, 'r') as infile:
        query_gene_alignments = [GenomicRangesAlignment.from_dict(dict_alignment)
                                 for dict_alignment in json.load(infile)]
    if os.path.exists(query_bg_filename):
        with open(query_bg_filename, 'r') as infile:
            query_bg = json.load(infile)
    else:
        query_bg = []
    try:
        background = fitter(query_bg)
        score_threshold = background.isf(pval_threshold)
        aligned = False
        query_orthologs = list()
        subject_orthologs = list()
        query_dropped = list()
        for alignment in query_gene_alignments:
            # Applies to FDR also to reduce memory consumption.
            alignment.filter_by_score(score_threshold)
            # Heavily optimized FDR_BH procedure.
            if fdr:
                pvals = background.sf([hsp.score
                                       for hsp in alignment.HSPs])
                total_hsps = len(alignment._all_HSPs)
                retain = pvals <= total_hsps * pval_threshold / ranking(pvals)
                del pvals
                alignment.filter_by_bool_array(retain)
            alignment.srange.name = alignment.qrange.name
            result = alignment.best_chain()
            try:
                record = result.to_bed12(mode='str')
                query_orthologs.append(record[0])
                subject_orthologs.append(record[1])
                aligned = True
            except ValueError:
                query_dropped = [alignment.qrange]
        if aligned:
            query_dropped = list()
        return query_orthologs, subject_orthologs, query_dropped
    except Exception as e:
        raise ExceptionLogger(e, query_gene_alignments[0].qrange)


build_orthologs_parser = argparse.ArgumentParser(description='Asses orthologous alignments based on chosen statistical strategy and build orthologs.',
                                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
build_orthologs_parser.add_argument('-alignments',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='Alignments dir produced with ortho2align get_alignments.')
build_orthologs_parser.add_argument('-background',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='background dir produced with ortho2align estimate_background.')
build_orthologs_parser.add_argument('-fitting',
                                    type=str,
                                    nargs='?',
                                    choices=['kde', 'hist'],
                                    default='kde',
                                    help='approach to fit background distribution (kde: KDE, hist: Histogram).')
build_orthologs_parser.add_argument('-threshold',
                                    type=float,
                                    nargs='?',
                                    default=0.05,
                                    help='p-value threshold to filter HSPs by score.')
build_orthologs_parser.add_argument('--fdr',
                                    action='store_true',
                                    help='use FDR correction for HSP scores.')
build_orthologs_parser.add_argument('-query_orthologs',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='output .bed12 filename for query orthologs.')
build_orthologs_parser.add_argument('-subject_orthologs',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='output .bed12 filename for subject orthologs.')
build_orthologs_parser.add_argument('-query_dropped',
                                    type=str,
                                    nargs='?',
                                    required=True,
                                    help='output .bed6 filename for dropped query ranges.')
build_orthologs_parser.add_argument('-query_exceptions',
                                    type=str,
                                    nargs='?',
                                    default='build_orthologs_exceptions.bed',
                                    help='output .bed6 filename for query ranges that caused exceptions while processing.')
build_orthologs_parser.add_argument('-cores',
                                    type=int,
                                    nargs='?',
                                    default=1,
                                    help='Number of cores to use for refinement multiprocessing.')
build_orthologs_parser.add_argument('-timeout',
                                    type=int,
                                    nargs='?',
                                    default=None,
                                    help='Time in seconds to terminate a single process of refinement of a single alignment.')
build_orthologs_parser.add_argument('--silent',
                                    action='store_true',
                                    help='silent CLI if included.')


def build_orthologs(alignments,
                    background,
                    fitting,
                    query_orthologs,
                    subject_orthologs,
                    query_dropped,
                    query_exceptions='build_orthologs_exceptions.bed',
                    threshold=0.05,
                    fdr=False,
                    cores=1,
                    timeout=None,
                    silent=False):
    from functools import partial
    from pathlib import Path

    from .fitting import HistogramFitter, KernelFitter
    from .genomicranges import BaseGenomicRangesList
    from .parallel import TimeoutProcessPool

    alignments_dir = Path(alignments)
    background_dir = Path(background)
    fitting_type = fitting
    query_output_filename = Path(query_orthologs)
    subject_output_filename = Path(subject_orthologs)
    query_dropped_filename = Path(query_dropped)
    query_exceptions_filename = Path(query_exceptions)
    pval_threshold = threshold

    cmd_hints = ['reading alignments...',
                 'getting background...',
                 'refining alignments and calling orthologs...',
                 'saving to files...',
                 'finished.']
    cmd_point = 0
    with tqdm(total=(len(cmd_hints) - 1),
              bar_format='{n_fmt}/{total_fmt} {elapsed}<{remaining} {postfix}',
              postfix=cmd_hints[cmd_point],
              disable=silent) as pbar:

        alignments_files = [os.path.join(alignments_dir, filename)
                            for filename in os.listdir(alignments_dir)
                            if filename.endswith('json')]
        bg_files = [os.path.join(background_dir, os.path.basename(alignments_file))
                    for alignments_file in alignments_files]

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()

        if fitting_type == 'hist':
            fitter = HistogramFitter
        elif fitting_type == 'kde':
            fitter = KernelFitter

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()

        build_partial = partial(_build,
                                fitter=fitter,
                                fdr=fdr,
                                pval_threshold=pval_threshold)
        data_for_building = zip(alignments_files, bg_files)

        with TimeoutProcessPool(cores, verbose=not silent) as p:
            try:
                orthologs, exceptions = p.starmap(build_partial,
                                                  data_for_building,
                                                  timeout)
            except Exception as e:
                raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')

        query_exception_ranges = list()
        if len(exceptions) > 0:
            tqdm.write(f'{len(exceptions)} exceptions occured:')
            for exception in exceptions:
                tqdm.write(str(exception))
                query_exception_ranges.append(exception.variable)

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()

        if len(orthologs) == 0:
            query_orthologs, subject_orthologs, query_dropped = [], [], []
        else:
            query_orthologs, subject_orthologs, query_dropped = zip(*orthologs)
        query_orthologs = [ortholog
                           for group in query_orthologs
                           for ortholog in group]
        subject_orthologs = [ortholog
                             for group in subject_orthologs
                             for ortholog in group]
        query_dropped = BaseGenomicRangesList([grange
                                               for group in query_dropped
                                               for grange in group])
        query_exception_list = BaseGenomicRangesList(query_exception_ranges)

        with open(query_output_filename, 'w') as outfile:
            for record in query_orthologs:
                outfile.write(record + "\n")
        with open(subject_output_filename, 'w') as outfile:
            for record in subject_orthologs:
                outfile.write(record + "\n")
        with open(query_dropped_filename, 'w') as outfile:
            query_dropped.to_bed6(outfile)
        with open(query_exceptions_filename, 'w') as outfile:
            query_exception_list.to_bed6(outfile)

        cmd_point += 1
        pbar.postfix = cmd_hints[cmd_point]
        pbar.update()


get_best_orthologs_parser = argparse.ArgumentParser(description='Select only one ortholog for each query gene based on provided variety of strategies.',
                                                    prog='ortho2align best_orthologs',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
get_best_orthologs_parser.add_argument('-query_orthologs',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='query orthologs bed12 file.')
get_best_orthologs_parser.add_argument('-subject_orthologs',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='subject orthologs bed12 file.')
get_best_orthologs_parser.add_argument('-value',
                                       type=str,
                                       nargs='?',
                                       choices=['total_length', 'block_count', 'block_length', 'weight'],
                                       required=True,
                                       help='which value of orthologs to use in case of multiple orthologs.')
get_best_orthologs_parser.add_argument('-function',
                                       type=str,
                                       nargs='?',
                                       choices=['max', 'min'],
                                       required=True,
                                       help='orthologs with which value to select in case of multiple orthologs.')
get_best_orthologs_parser.add_argument('-outfile_query',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='output filename for query orthologs.')
get_best_orthologs_parser.add_argument('-outfile_subject',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='output filename for subject orthologs.')
get_best_orthologs_parser.add_argument('-outfile_map',
                                       type=str,
                                       nargs='?',
                                       required=True,
                                       help='output json filename for mapping of query and subject ortholog names.')


def get_best_orthologs(query_orthologs,
                       subject_orthologs,
                       value,
                       function,
                       outfile_query,
                       outfile_subject,
                       outfile_map):
    import json
    from collections import defaultdict
    from pathlib import Path

    query_filename = Path(query_orthologs)
    subject_filename = Path(subject_orthologs)
    value_name = value
    func_name = function
    query_output_filename = Path(outfile_query)
    subject_output_filename = Path(outfile_subject)
    outfile_map_filename = Path(outfile_map)

    extract = {'total_length': lambda line: (int(line.strip().split('\t')[2])
                                             - int(line.strip().split('\t')[1])),
               'block_count': lambda line: int(line.strip().split('\t')[9]),
               'block_length': lambda line: sum([int(i)
                                                 for i in line.strip().split('\t')[10].split(',')]),
               'weight': lambda line: float(line.strip().split('\t')[4]),
               'name': lambda line: line.strip().split('\t')[3]}

    func = {'min': min, 'max': max}

    name_pairs = list()

    with open(query_filename, 'r') as query_in, \
         open(subject_filename, 'r') as subject_in, \
         open(query_output_filename, 'w') as query_out, \
         open(subject_output_filename, 'w') as subject_out:
        purge = list()
        for query_line, subject_line in zip(query_in, subject_in):
            if len(purge) == 0:
                purge.append((query_line, subject_line))
            elif extract['name'](purge[-1][0]) == extract['name'](query_line):
                purge.append((query_line, subject_line))
            else:
                best_ortholog = func[func_name](purge,
                                                key=lambda i: extract[value_name](i[0]))
                query_out.write(best_ortholog[0])
                subject_out.write(best_ortholog[1])
                name_pairs.append([extract['name'](best_ortholog[0]),
                                   extract['name'](best_ortholog[1])])
                purge = list()
                purge.append((query_line, subject_line))

        if len(purge) > 0:
            best_ortholog = func[func_name](purge,
                                            key=lambda i: extract[value_name](i[0]))
            query_out.write(best_ortholog[0])
            subject_out.write(best_ortholog[1])
            name_pairs.append([extract['name'](best_ortholog[0]),
                               extract['name'](best_ortholog[1])])

    the_map = defaultdict(list)
    for query_name, subject_name in name_pairs:
        the_map[query_name].append(subject_name)
    with open(outfile_map_filename, 'w') as outfile:
        json.dump(the_map, outfile)


benchmark_orthologs_parser = argparse.ArgumentParser(description='Compare found orthologs against real orthologs and calculate several performance metrics.',
                                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
benchmark_orthologs_parser.add_argument('-query_genes',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='query genomic ranges.')
benchmark_orthologs_parser.add_argument('-query_name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the query genes annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
benchmark_orthologs_parser.add_argument('-found_query',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='found orthologs of query genomic ranges in query genome.')
benchmark_orthologs_parser.add_argument('-found_query_name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the found query orthologs annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
benchmark_orthologs_parser.add_argument('-found_subject',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='found orthologs of query genomic ranges in subject genome.')
benchmark_orthologs_parser.add_argument('-found_subject_name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the found subject orthologs annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
benchmark_orthologs_parser.add_argument('-found_query_map',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='json map linking query genes names and names of corresponding found query orthologs.')
benchmark_orthologs_parser.add_argument('-found_subject_map',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='json map linking query genes names and names of corresponding found subject orthologs.')
benchmark_orthologs_parser.add_argument('-found_query_subject_map',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='json map linking query orthologs names and corresponding subject orthologs names.')
benchmark_orthologs_parser.add_argument('-real_subject',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='real orthologs of query genomic ranges in subject genome.')
benchmark_orthologs_parser.add_argument('-real_subject_name_regex',
                                        type=str,
                                        nargs='?',
                                        default=None,
                                        help='Regular expression for extracting gene names from the real subject orthologs annotation (.gff and .gtf only). '
                                             'Must contain one catching group.')
benchmark_orthologs_parser.add_argument('-real_map',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='json map linking query genes names and names of corresponding real orthologs.')
benchmark_orthologs_parser.add_argument('-outfile',
                                        type=str,
                                        nargs='?',
                                        required=True,
                                        help='json output filename.')
benchmark_orthologs_parser.add_argument('-tp_mode',
                                        type=str,
                                        nargs='?',
                                        choices=['all', 'single'],
                                        default='all',
                                        help='how to calculate true positives')


def benchmark_orthologs(query_genes,
                        query_name_regex,
                        found_query,
                        found_query_name_regex,
                        found_subject,
                        found_subject_name_regex,
                        found_query_map,
                        found_subject_map,
                        found_query_subject_map,
                        real_subject,
                        real_subject_name_regex,
                        real_map,
                        outfile,
                        tp_mode):
    import json

    from .benchmark import calc_ortholog_metrics, trace_orthologs
    from .parsing import parse_annotation

    query_genes_filename = query_genes
    found_query_filename = found_query
    found_subject_filename = found_subject
    found_query_map_filename = found_query_map
    found_subject_map_filename = found_subject_map
    found_query_subject_map_filename = found_query_subject_map
    real_subject_filename = real_subject
    real_map_filename = real_map
    out_filename = outfile

    with open(query_genes_filename, 'r') as infile:
        query_genes = parse_annotation(infile,
                                       name_regex=query_name_regex)

    with open(found_query_filename, 'r') as infile:
        found_query_genes = parse_annotation(infile,
                                             name_regex=found_query_name_regex)

    with open(found_subject_filename, 'r') as infile:
        found_subject_genes = parse_annotation(infile,
                                               name_regex=found_subject_name_regex)

    with open(real_subject_filename, 'r') as infile:
        real_subject_genes = parse_annotation(infile,
                                              name_regex=real_subject_name_regex)

    with open(found_query_map_filename, 'r') as infile:
        found_query_map = json.load(infile)

    with open(found_subject_map_filename, 'r') as infile:
        found_subject_map = json.load(infile)

    with open(found_query_subject_map_filename, 'r') as infile:
        found_query_subject_map = json.load(infile)

    with open(real_map_filename, 'r') as infile:
        real_map = json.load(infile)

    query_genes.relation_mapping(found_query_genes, found_query_map, 'found_query')
    query_genes.relation_mapping(found_subject_genes, found_subject_map, 'found_subject')
    query_genes.relation_mapping(real_subject_genes, real_map, 'real')

    for query_gene in query_genes:
        query_gene.relations['found_query'].relation_mapping(query_gene.relations['found_subject'],
                                                             found_query_subject_map,
                                                             'ortho_link')

    trace_orthologs(query_genes)
    metrics = calc_ortholog_metrics(query_genes, tp_mode=tp_mode)

    with open(out_filename, 'w') as outfile:
        json.dump(metrics, outfile)


def run_pipeline(config):
    # validate the config
    # combine functions above with data from the config
    pass