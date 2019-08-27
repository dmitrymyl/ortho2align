import pandas as pd
from pathlib import Path
import os


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
