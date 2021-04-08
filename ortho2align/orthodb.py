import json
import os
from collections import defaultdict
from pathlib import Path

import pandas as pd
import tqdm


def load_odb_file(file_suffix,
                  odb_path,
                  odb_prefix,
                  colnames,
                  dtype,
                  **kwargs):
    """Loads OrthoDB table into pandas.DataFrame.

    The loaded file is `odb_path + odb_prefix + file_suffix`.

    Args:
        file_suffix (str): suffix of OrthoDB table to be loaded.
        odb_path (str, Path): path to OrthoDB folder.
        odb_prefix (str): prefix of OrthoDB files.
        colnames (list): list of column names of OrthoDB table.
        dtype (dict): dictionary of column dtypes as {'colname': 'coldtype'}.
        **kwargs: any other arguments to be passed to `pandas.read_csv()`.

    Returns:
        pandas.DataFrame instance with loaded OrthoDB table.
        Columns are named as in `colnames`, column dtypes are
        as in `dtypes`.

    Raises:
        ValueError in case file doesn't exist.
    """

    if not isinstance(odb_path, Path):
        odb_path = Path(odb_path)

    path = odb_path.joinpath(odb_prefix + file_suffix)
    if path.exists():
        return pd.read_csv(str(path),
                           header=None,
                           sep='\t',
                           names=colnames,
                           dtype=dtype,
                           **kwargs)
    else:
        raise ValueError(f"The path {path} doesn't exist.")


def load_table(path, colnames, dtype, **kwargs):
    """Loads tab-separated table into memory.

    Args:
        path (str, Path): path to table.
        colnames (list): list of column names of OrthoDB table.
        dtype (dict): dictionary of column dtypes as {'colname': 'coldtype'}.
        **kwargs: any other arguments to be passed to `pandas.read_csv()`.

    Returns:
        pandas.DataFrame instance with loaded table.
        Columns are named as in `colnames`, column dtypes are
        as in `dtypes`.
    """
    return pd.read_csv(str(path),
                       header=None,
                       sep='\t',
                       names=colnames,
                       dtype=dtype,
                       **kwargs)


def query_table(infile, outfile, column, term, column_list, sep):
    """Queries table column by term.

    Finds all lines in infile, that contain
    `term` in column named `column` and writes
    that lines to outfile with excluded column
    that was queried on.

    Args:
        infile (file-like): input file.
        outfile (file-like): output file.
        column (str): column name to query on.
        term (str): a term to find in column.
        column_list (list): list of column names.
        sep (str): field separator.

    Returns:
        None.
    """
    lookup_index = column_list.index(column)
    for line in infile:
        lst = line.strip().split(sep)
        if lst[lookup_index] == term:
            lst.pop(lookup_index)
            outfile.write(sep.join(lst) + '\n')


def query_cached_odb_file(file_suffix,
                          odb_path,
                          odb_prefix,
                          column,
                          term,
                          column_list,
                          cache_path):
    """Queries cached OrthoDB file.

    In case of the first query of the given
    file in given column on given term
    it is queried and written to cache
    directory.

    Args:
        file_suffix (str): suffix of OrthoDB table to be loaded.
        odb_path (str, Path): path to OrthoDB folder.
        odb_prefix (str): prefix of OrthoDB files.
        column (str): column name to query on.
        term (str): a term to find in column.
        column_list (list): list of column names.
        cache_path (str, Path): a path to cache directory.

    Returns:
        (Path) a path to queried file in the following scheme:
        cache_path + "/" + odb_prefix + odb_suffix +
        "_" + column + "_" + term
    """
    if not isinstance(odb_path, Path):
        odb_path = Path(odb_path)
    if not isinstance(cache_path, Path):
        cache_path = Path(cache_path)

    cache_filename = cache_path.joinpath("".join([odb_prefix,
                                                  file_suffix,
                                                  "_",
                                                  column,
                                                  "_",
                                                  term]))
    if not cache_filename.exists():
        infilename = odb_path.joinpath(odb_prefix + file_suffix)
        with open(infilename, 'r') as infile:
            if not cache_path.exists():
                os.mkdir(cache_path)
            with open(cache_filename, 'w') as outfile:
                query_table(infile, outfile, column, term, column_list, '\t')
    return cache_filename


def filter_table(infile,
                 outfile,
                 column,
                 terms,
                 column_list,
                 sep):
    """Filters table by terms in given column.

    Leaves only those lines that contain one of
    the terms in specified column.

    Args:
        infile (file-like): input file.
        outfile (file-like): output file.
        column (str): column name to query on.
        terms (set): a set of terms to find in column.
        column_list (list): list of column names.
        sep (str): field separator.

    Returns:
        None.
    """
    lookup_index = column_list.index(column)
    for line in infile:
        lst = line.strip().split(sep)
        if lst[lookup_index] in terms:
            outfile.write(line)


def filter_odb_file(file_suffix,
                    odb_path,
                    odb_prefix,
                    column,
                    terms,
                    column_list,
                    cache_path,
                    code):
    """Filters OrthoDB table by terms in given column.

    Leaves only those lines that contain one of
    the terms in specified column.

    Args:
        file_suffix (str): suffix of OrthoDB table to be filtered.
        odb_path (str, Path): path to OrthoDB folder.
        odb_prefix (str): prefix of OrthoDB files.
        column (str): column name to filter on.
        terms (set): a set of terms to find in column.
        column_list (list): list of column names.
        cache_path (str, Path): a path to cache directory.
        code (str): a code to write at the end of the filtered
            file.

    Returns:
        (Path) a path to queried file in the following scheme:
        cache_path + "/" + odb_prefix + odb_suffix +
        "_" + column + "-" + code
    """

    if not isinstance(odb_path, Path):
        odb_path = Path(odb_path)
    if not isinstance(cache_path, Path):
        cache_path = Path(cache_path)

    cache_filename = cache_path.joinpath("".join([odb_prefix,
                                                  file_suffix,
                                                  "_",
                                                  column,
                                                  "-",
                                                  code]))
    if not cache_path.exists():
        os.mkdir(cache_path)
    infilename = odb_path.joinpath(odb_prefix + file_suffix)
    with open(infilename, 'r') as infile:
        with open(cache_filename, 'w') as outfile:
            filter_table(infile, outfile, column, terms, column_list, '\t')
    return cache_filename


class FileOperationWrapper:
    """Provides context manager and dictionary-
    based file handler for split_table function.

    Attributes:
        file_dict (dict): handler of opened files.
        prefix (str, Path): file prefix for each of
            the opened files.
    """

    def __init__(self, prefix):
        """Initializes FileOperationWrapper object.

        Args:
            prefix (str, Path): file prefix.

        Returns:
            None.
        """
        self.file_dict = dict()
        self.prefix = prefix

    def __setitem__(self, key, item):
        """Item setter.

        Args:
            key (str): a key to file instance.
            item (file-like): file instance.

        Returns:
            None.
        """
        self.file_dict[key] = item

    def __getitem__(self, key):
        """Item getter.

        In case file associated with the key
        doesn't exist in `self.file_dict`, it
        is opened for writing and added to
        `self.file_dict`.

        Args:
            key (str): a key to file instance.

        Returns:
            file instance opened at the location
            with following scheme:
                self.prefix + "_" + key
        """
        if key not in self.file_dict.keys():
            new_file_name = str(self.prefix) + "_" + key
            self.file_dict[key] = open(new_file_name, 'w')
        return self.file_dict[key]

    def __enter__(self):
        return self

    def __exit__(self, ex_type, ex_value, ex_traceback):
        for file in self.file_dict.values():
            file.close()


def split_table(infile, outprefix, column, column_list, sep):
    """Splits table on many files based on terms in specific column.

    For each unique term in column a file is created
    and lines containing specific term in the column
    are written to corresponding files. The column that
    was filtered on is omitted. New files are named
    in the following scheme:
    outprefix + "_" + term.

    Args:
        infile (file-like): input file.
        outprefix (str, Path): output files prefix.
        column (str): column name to filter on.
        column_list (list): list of column names.
        sep (str): field separator.

    Returns:
        None.
    """
    lookup_index = column_list.index(column)

    with FileOperationWrapper(outprefix) as file_handler:
        for line in infile:
            lst = line.strip().split(sep)
            split_term = lst[lookup_index]
            lst.pop(lookup_index)
            file_handler[split_term].write(sep.join(lst) + "\n")


def split_odb_file(file_suffix,
                   odb_path,
                   odb_prefix,
                   column,
                   column_list,
                   cache_path):
    """Splits OrthoDB table on many files based on terms in specific column.

    For each unique term in column a file is created
    and lines containing specific term in the column
    are written to corresponding files. The column that
    was filtered on is omitted. New files are named
    in the following scheme:
    cache_path "/" + odb_prefix + file_suffix + "_" +
    column + "_" + term.

    Args:
        file_suffix (str): suffix of OrthoDB table to be splitted.
        odb_path (str, Path): path to OrthoDB folder.
        odb_prefix (str): prefix of OrthoDB files.
        column (str): column name to split on.
        column_list (list): list of column names.
        cache_path (str, Path): a path to cache directory.

    Returns:
        None.
    """
    if not isinstance(odb_path, Path):
        odb_path = Path(odb_path)
    if not isinstance(cache_path, Path):
        cache_path = Path(cache_path)

    infilename = odb_path.joinpath(odb_prefix + file_suffix)
    if not cache_path.exists():
        os.mkdir(cache_path)
    outprefix = cache_path.joinpath("".join([odb_prefix,
                                             file_suffix,
                                             "_",
                                             column]))
    with open(infilename, 'r') as infile:
        split_table(infile, outprefix, column, column_list, '\t')


def cache_orthodb_xrefs(orthodb_path, cache_path, orthodb_prefix, xref_suffix):
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


def get_orthodb_map(query_genes, output_json_file, query_db, subject_db,
                    level_taxid, subject_taxids, orthodb_path, cache_path,
                    orthodb_prefix, silent, tmp_path):
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
