import pandas as pd
from pathlib import Path
import os


def load_odb_file(file_suffix,
                  odb_path,
                  odb_prefix,
                  colnames,
                  dtype,
                  **kwargs):

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
    return pd.read_csv(str(path),
                       header=None,
                       sep='\t',
                       names=colnames,
                       dtype=dtype,
                       **kwargs)


def query_table(infile, outfile, column, term, column_list, sep):
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
            file instance.
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
    lookup_index = column_list.index(column)

    with FileOperationWrapper(outprefix) as file_handler:
        for line in infile:
            lst = line.strip().split(sep)
            split_term = lst[lookup_index]
            lst.pop(lookup_index)
            file_handler[split_term].write(sep.join(lst) + "\n")


def split_odb_table(file_suffix,
                    odb_path,
                    odb_prefix,
                    column,
                    column_list,
                    cache_path):

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
