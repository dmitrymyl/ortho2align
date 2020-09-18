"""Chromosome coordinate systems.
0-based: chromosome start coordinate is 0.
1-based: chromosome start coordinate is 1.
inclusive: character under end coordinate
    is included in the sequence.
exclusive: character under end coordinate
    is not included in the sequence.

Here we use 0-based exclusive system for
easier arithmetics.

Translation from other systems:
0-based inclusive: start += 0, end += 1
1-based inclusive: start += -1, end += 0
1-based exclusive: start += -1, end += -1
"""
import re
from itertools import chain
from pathlib import Path
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import namedtuple, defaultdict
from sortedcontainers import SortedKeyList
from tqdm import tqdm
from .alignment import HSPVertex, Alignment, AlignmentChain, nxor_strands, numberize


class GenomicException(Exception):
    """Basic exception in genomicranges module."""
    pass


class ChromosomeNotFoundError(GenomicException):
    """Occurs when chromosome was not found in given genome."""

    def __init__(self, chrom, genome, list_of_chroms):
        self.chrom = chrom
        self.genome = genome
        self.list_of_chroms = list_of_chroms
        super().__init__(f"Chromosome {self.chrom} is not "
                         f"in {self.genome} among these "
                         f"chromosomes: {', '.join(list_of_chroms)}")


class InconsistentChromosomesError(GenomicException):
    """Occurs when two GenomicRange instances have different chromosomes."""

    def __init__(self, one_chrom, other_chrom):
        self.one_chrom = one_chrom
        self.other_chrom = other_chrom
        super().__init__(f"Inconsistent chromosomes of "
                         f"comparing genomic ranges: "
                         f"{self.one_chrom} and {self.other_chrom}.")


class InconsistentGenomesError(GenomicException):
    """Occurs when two GenomicRange instances have different genomes."""

    def __init__(self, one_genome, other_genome):
        self.one_genome = one_genome
        self.other_genome = other_genome
        super().__init__(f"Inconsistent genomes of "
                         f"comparing genomic ranges: "
                         f"{self.one_genome} and {self.other_genome}.")


class EmptyGenomicRangesListError(GenomicException):
    """Occurs when GenomicRangesList is empty."""

    def __init__(self, range_list_instance):
        self.range_list = range_list_instance
        super().__init__(f"Empty GenomicRangesList instance: "
                         f"{self.range_list}")


class GenomicCoordinatesError(GenomicException):
    """Occurs when coordinates of GenomicRange instance are out of the chromosome."""

    def __init__(self, grange, chromosome):
        self.grange = grange
        self.chromosome = chromosome
        super().__init__(f"Coordinates of GenomicRange {self.grange} "
                         f"are out of chromosome with size "
                         f"{self.chromosome.size}")


class FilePathNotSpecifiedError(GenomicException):
    """Occurs when location of sequence file is not specified, but is accessed."""

    def __init__(self, instance):
        self.instance = instance
        super().__init__(f"Sequence file location for {self.instance} "
                         f"is not specified.")


class SequenceFileNotFoundError(GenomicException):
    """Occurs when provided sequence file path is invalid."""

    def __init__(self, instance, file_path):
        self.instance = instance
        self.file_path = file_path
        super().__init__(f"Sequence file {self.file_path} "
                         f"for {self.instance} doesn't exist.")


from_line = 'AaTtGgCc'
to_line = 'TtAaCcGg'
trans_map = str.maketrans(from_line, to_line)


def fasta_getter(infile, outfile, total=None, linewidth=60,
                 reverse=False, complement=False):
    """Reformats fasta file.

    Allows one to reformat existing sequence in
    fasta format to the specified linewidth. The pointer in the infile
    must be put at the beginning of the desired sequence.

    Args:
        infile (fileobj): input fasta file with cursor set
            at the start of sequence fragment.
        outfile (fileobj): output fasta file.
        total (int): number of symbols to read (default: None;
            means to read all sequence up to EOF or next record).
        linewidth (int): linewidth of new sequence
            (default: 60).
        reverse (bool): if True, takes sequence backwards
            (default: False).
        complement (bool): if True, returns complement
            nucleotide sequence (default: False).
    Returns:
        None
    """
    # Getting the record.
    record = infile.read(total).replace("\n", "")
    # Processing the recod.
    if complement:
        record = record.translate(trans_map)
    if reverse:
        record = record[::-1]
    # writing the record
    record_to_write = "\n".join([record[i: i + linewidth]
                                 for i in range(0, len(record), linewidth)] + [""])
    outfile.write(record_to_write)


class SequencePath:
    """Wrapper aroudn pathlib.Path with checking of None.

    Attributes:
        path (None or Path): filesystem path.
    """
    __slots__ = ('path',)

    def __init__(self, path=None):
        """Initializes SequencePath instance.

        Args:
            path (str, bytes, Path): filesystem path
                (default: None).

        Returns:
            None
        """
        self.path = None if path is None else Path(path)

    def __str__(self):
        """Returns str representation."""
        return str(self.path)

    def __repr__(self):
        """Returns representation."""
        return f"SequencePath('{self.path}')"

    def __fspath__(self):
        """Magic method to allow use with built-in `open()`."""
        return str(self.path)

    def __eq__(self, other):
        """Equality magic method."""
        return self.path == other.path

    def __truediv__(self, key):
        if self.path is None:
            return None
        else:
            return self.__class__(self.path / key)

    def __rtruediv__(self, key):
        if self.path is None:
            return None
        else:
            return self.__class(key / self.path)

    def to_dict(self):
        """Returns dict representation."""
        if self.path is None:
            return {'path': None}
        else:
            return {'path': str(self.path)}

    @classmethod
    def from_dict(cls, dict_):
        return cls(dict_.get('path'))

    def check_correct(self, foreign=None):
        """Checks whether the path exists or it is correct.

        Checks whether the path is not None, whether it
        exists and whether it is a Path instance. If None
        or doesn't exist, reports. If not a Path instance,
        casts to the Path.

        Args:
            foreign (object): foreign class instance
                that contains SequencePath instance
                (default: None).

        Returns:
            None

        Raises:
            FileNotSpecifiedError in case self.path is None.
            SequenceFileNotFoundError in case self.path
                points to the path that doesn't exist.
        """
        if self.path is None:
            raise FilePathNotSpecifiedError(foreign)
        if not isinstance(self.path, Path):
            self.path = Path(self.path)
        if not self.path.exists():
            raise SequenceFileNotFoundError(foreign, self.path)

    def exists(self):
        """Checks whether self.path exists."""
        if self.path is None:
            return False
        return self.path.exists()

    def open(self, mode='r'):
        """Opens corresponding file.

        Args:
            mode (str): mode to open file
                (default: 'r').

        Returns:
            (fileobj) opened file instance.
        """
        return open(self.path, mode=mode)


class BaseGenomicRange:
    """
    Base class for any genomic ranges.

    Genomic coordinates are 0-based exclusive:
    * chromosome start coordinate is 0;
    * character under genomic range coordinate
        is not included into the sequence.

    Uses `__slots__` instead of `__dict__` for fast
    attribute access.

    Implemented slots:
        chrom
        start
        end
        strand
        _name

    Attributes:
        chrom (str, int): name of the chromosome.
        start (int): start of the genomic range.
        end (int): end of the genomic range.
        strand (str): strand of the genomic range.
            Must be one of `"+"`, `"-"` or `"."` (inessential
            strandness) (default: `"."`).
        _name (str, int): "private" variable for name handling (default: `None`).
        sequence_header (str): a property representing genomic range
            coordinates in POSIX-compatible way of file names.
        name (str): property representing the name of the genomic range.
            if `_name` is `None`, then falls for `sequence_header`, otherwise
            it is `_name`.

    Basic methods:
        __init__
        __len__
        __repr__
        __str__
        __eq__

    Properties:
        name
        sequence_header
        init_args

    Dictionary interoperability:
        to_dict
        from_dict

    Genomic methods:
        distance
        intersect
        merge
        flank
        find_neighbours

    Benchmarking methods:
        calc_JI
        calc_OC
        calc_fraction

    When subclassing, following methods should be
    reimplemented with super() call:
        __repr__
        __str__
        __eq__
        to_dict
        from_dict
        distance
        intersect
        merge
        flank
        find_neighbours
    """
    __slots__ = ('chrom', 'start', 'end', 'strand', '_name')

    def __init__(self, chrom, start, end, strand='.', name=None):
        """Initializes BaseGenomicRange instance.

        Args:
            chrom (str, int): name of the chromosome.
            start (int): start of the genomic range.
            end (int): end of the genomic range.
            strand (str): strand of the genomic range,
                one of `"+"`, `"-"`, `"."` (default: `"."`).
            name (str, int): name of the genomic range
                (default: `None`).

        Returns:
            None
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self._name = name

    def __len__(self):
        """Returns the length of the genomic range.

        Returns:
            int: length of the genomic range.
        """
        return self.end - self.start

    @property
    def init_args(self):
        """Returns tuple of __init__ arguments."""
        return 'chrom', 'start', 'end', 'strand', 'name'

    def __repr__(self):
        """Returns code representation of the genomic range."""
        attrs = ", ".join((repr(getattr(self, attr))
                           for attr in self.init_args))
        return f"BaseGenomicRange({attrs})"

    def __str__(self):
        """Returns tab-separated string representation of the genomic range."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.strand}"

    def __eq__(self, other):
        """Checks whether two genomic ranges are equal.

        Args:
            other (BaseGenomicRange): the other genomic range.

        Returns:
            bool: whether two genomic ranges have equal attributes.
        """
        return all((getattr(self, attr) == getattr(other, attr)
                    for attr in self.init_args))

    def to_dict(self):
        """Returns dictionary representation of the genomic range.

        Might be useful for saving to json file format.

        Returns:
            dict: dictionary representation of the genomic range.
        """
        return {'chrom': self.chrom,
                'start': self.start,
                'end': self.end,
                'strand': self.strand,
                'name': self._name}

    @classmethod
    def from_dict(cls, dict_):
        """Gets genomic range from its dictionary representation.

        Args:
            dict_ (dict): a dictionary generated with
                `BaseGenomicRange.to_dict`.

        Returns:
            BaseGenomicRange: restored genomic range.
        """
        parser_dict = {'chrom': lambda i: i.get('chrom'),
                       'start': lambda i: i.get('start'),
                       'end': lambda i: i.get('end'),
                       'strand': lambda i: i.get('strand'),
                       'name': lambda i: i.get('name')}
        return cls(**{attr: func(dict_)
                      for attr, func in parser_dict.items()})

    @property
    def sequence_header(self):
        """Returns sequence header in POSIX-compatible file name system."""
        return f"{self.chrom}:{self.start}-{self.end}{self.strand}"

    @property
    def name(self):
        """Returns `_name` or `sequence_header` if former is `None`."""
        if self._name is None:
            return self.sequence_header
        return self._name

    @name.setter
    def name(self, value):
        """Sets new `name` value through updating `_name`."""
        self._name = value

    def distance(self, other):
        """Calculates distance between two genomic ranges.

        In case chromosomes are consistent, returns distance
        between two genomic ranges in strand-nonspecific way.

        Args:
            other (BaseGenomicRange or childs): the other genomic range
                to calculate distance to `self`.

        Raises:
            InconsistentChromosomesError: an error in case of inconsistent chromosomes.
        """
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        if self.end < other.start:
            return other.start - self.end
        if other.end < self.start:
            return other.end - self.start
        return 0

    def intersect(self, other):
        """Intersects two genomic ranges.

        In case distance between two genomic ranges equals 0,
        returns intersection of them in strand-nonspecific way.
        In the other case returns `None`.

        Args:
            other (BaseGenomicRange): the other genomic range
                to be intersected with `self`.

        Returns:
            BaseGenomicRange or None: an intersection of two genomic ranges
            with `"."` strand.

        Raises:
            InconsistentChromosomesError: an error in case of inconsistent chromosomes.
        """
        if self.distance(other) == 0:
            if self.start < other.start < self.end:
                start = other.start
            else:
                start = self.start
            if self.start < other.end < self.end:
                end = other.end
            else:
                end = self.end
            used_args = {'chrom', 'start', 'end', 'strand', 'name'}
            kwargs = {attr: getattr(self, attr)
                      for attr in set(self.init_args) - used_args}
            return self.__class__(chrom=self.chrom,
                                  start=start,
                                  end=end,
                                  strand='.',
                                  **kwargs)
        return None

    def calc_JI(self, other):
        """Calculates Jaccard index between two genomic ranges.

        Jaccard index (JI) is the ratio of the length of the intersection
        to the length of the union of two genomic ranges. In case division
        by zero encountered, returns zero.

        Args:
            other (BaseGenomicRange): the other genomic range.

        Returns:
            float: JI value.

        Raises:
            InconsistentChromosomesError: an error in case of inconsistent chromosomes.
        """
        try:
            intersection = self.intersect(other)
            if intersection is None:
                return 0
            return len(intersection) / (len(self) + len(other) - len(intersection))
        except GenomicException:
            return 0
        except ZeroDivisionError:
            return 0

    def calc_OC(self, other):
        """Calculates overlap coefficient between two genomic ranges.

        Overlap coefficient (OC) is the ratio of the length of the
        intersection to the length of the shortest grange of the two
        presented. In case division by zero encountered, returns zero.

        Args:
            other (BaseGenomicRange): the other genomic range.

        Returns:
            float: OC value.

        Raises:
            InconsistentChromosomesError: an error in case of inconsistent chromosomes.
        """
        try:
            intersection = self.intersect(other)
            if intersection is None:
                return 0
            return len(intersection) / min(len(self), len(other))
        except GenomicException:
            return 0
        except ZeroDivisionError:
            return 0

    def calc_fraction(self, other):
        """Calculates fraction of `self` overlapped with `other`.

        Returns the ratio of the length of the intersection
        of two genomic ranges to the length of `self`. In case division 
        by zero encountered, returns zero.

        Args:
            other (BaseGenomicRange): the other genomic range.

        Returns:
            float: fraction of `self` overlapped with `other`.

        Raises:
            InconsistentChromosomesError: an error in case of inconsistent chromosomes.
        """
        try:
            intersection = self.intersect(other)
            if intersection is None:
                return 0
            return len(intersection) / len(self)
        except GenomicException:
            return 0
        except ZeroDivisionError:
            return 0

    def merge(self, other):
        """Merges two genomic ranges together.

        Returns a new genomic range with minimal
        start coordinate and maximal end coordinate.
        Strand resets to `"."`.

        Returns:
            BaseGenomicRange: merge result.

        Raises:
            InconsistenChromosomesError: an error in case of inconsistent chromosomes.
        """
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        used_args = {'chrom', 'start', 'end', 'strand', 'name'}
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(chrom=self.chrom,
                              start=start,
                              end=end,
                              strand=".",
                              **kwargs)

    def flank(self, flank_distance=0):
        """Flanks genomic range with the given distance.

        Start coordinate is decreased and end
        coordinate is increased by the flank distance.
        Do not check if the start coordinate is negative or
        the end coordinate exceeds chromosome size.

        Args:
            flank_distance (int): distance to flank
                genomic range.

        Returns:
            BaseGenomicRange: flanked genomic range.

        Raises:
            ValueError: an error in case `flank_distance` value
                is negative and greater than the half-length
                of the genomic range.
        """
        if - flank_distance * 2 >= len(self):
            raise ValueError(f"Flank distance value {flank_distance} "
                             "is negative and its modulus is equal to "
                             "or greater than the half-length of the "
                             "genomic range.")
        used_args = {'chrom', 'start', 'end', 'strand', 'name'}
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(chrom=self.chrom,
                              start=self.start - flank_distance,
                              end=self.end + flank_distance,
                              strand=".",
                              **kwargs)

    def find_neighbours(self, other_list, distance=0):
        """Finds neighbours of self in other list at given distance.

        Args:
            other_list (BaseGenomicRangesList, GenomicRangesList): other list
                to find neighbours in.
            distance (int, float, default 0): distance at which two
                genomic ranges are considered neighbours.

        Returns:
            ``type(other_list)``: genomic ranges list with all
            neighbours of self in other_list at given distance.
        """
        neighbours = list()
        index = other_list.bisect_left(self)
        left_index = index - 1
        right_index = index
        try:
            while True and left_index >= 0:
                if abs(self.distance(other_list[left_index])) <= distance:
                    neighbours.append(other_list[left_index])
                left_index -= 1
        except GenomicException:
            pass
        except IndexError:
            pass
        try:
            while abs(self.distance(other_list[right_index])) <= distance:
                neighbours.append(other_list[right_index])
                right_index += 1
        except GenomicException:
            pass
        except IndexError:
            pass
        used_args = {'collection', }
        kwargs = {attr: getattr(other_list, attr)
                  for attr in set(other_list.init_args) - used_args}
        return other_list.__class__(neighbours, **kwargs)


class GenomicRange(BaseGenomicRange):
    """Represents genomic range for real applications.

    Inherits from `BaseGenomicRange` and extends its
    slots with *genome*, *_sequence_file_path* and
    *relations*.

    Attributes:
        chrom (str, int): Name of the chromosome.
        start (int): Start of the genomic range.
        end (int): End of the genomic range.
        strand (str): Strand of the genomic range,
            one of *"+"*, *"-"*, *"."*.
        name (str, int): Name of the genomic range.
        genome (str): Corresponding genome filename.
        sequence_file_path (SequencePath, str): Location
            of the genomic range sequence file.
        relations (dict): A dictionary with related
            `BaseGenomicRangesList` or `GenomicRangesList` instances.
        init_args (tuple): Returns a tuple of essential *__init__* arguments.
            This property is used for a convenient creation of
            new instances of the class in inherited methods.
    """

    __slots__ = ('genome', '_sequence_file_path', 'relations')

    def __init__(self, chrom, start, end, strand='.', name=None, genome=None,
                 sequence_file_path=None, relations=None, **kwargs):
        """Initializes GenomicRange instance.

        In practice, *genome*, *sequence_file_path* and *relations*
        should be set to ``None``. *genome* is used when parsing
        annotation files, *sequence_file_path* is used when
        extracting sequences via `FastaSeqFile` with `GenomicRangesList`,
        *relations* are used in several applications. *kwargs* is just
        a convenient placeholder.

        Args:
            chrom (str, int): Name of the chromosome.
            start (int): Start of the genomic range.
            end (int): End of the genomic range.
            strand (str, default "."): Strand of the genomic range,
                one of *"+"*, *"-"*, *"."*.
            name (str, int, default None): Name of the genomic range.
            genome (str, default None): Corresponding genome filename.
            sequence_file_path (SequencePath, str, default None): Location
                of the genomic range sequence file.
            relations (dict, default None): A dictionary with related
                BaseGenomicRangesList or GenomicRangesList instances.
            kwargs (dict): Other arguments.


        Returns:
            None
        """
        super().__init__(chrom=chrom,
                         start=start,
                         end=end,
                         strand=strand,
                         name=name)

        self.genome = genome
        self._sequence_file_path = (sequence_file_path
                                    if isinstance(sequence_file_path,
                                                  SequencePath)
                                    else SequencePath(sequence_file_path))
        self.relations = dict() if relations is None else relations

    @property
    def init_args(self):
        return 'chrom', 'start', 'end', 'strand', 'name', 'genome'

    def __repr__(self):
        """Returns code representation of the genomic range."""
        attrs = ('chrom',
                 'start',
                 'end',
                 'strand',
                 'name',
                 'genome',
                 'sequence_file_path',
                 'relations')
        return f"GenomicRange({', '.join([repr(getattr(self, attr)) for attr in attrs])})"

    def __str__(self):
        """Returns tab-separated table representation of the genomic range."""
        attrs = ('chrom',
                 'start',
                 'end',
                 'name',
                 'strand',
                 'genome')
        return "\t".join(str(getattr(self, attr))
                         for attr in attrs)

    def __eq__(self, other):
        """Checks for equality of two genomic ranges."""
        attrs = ('chrom',
                 'start',
                 'end',
                 'strand',
                 'name',
                 'genome',
                 'sequence_file_path',
                 'relations')
        return all(getattr(self, attr) == getattr(other, attr)
                   for attr in attrs)

    def to_dict(self):
        """Returns dictionary representation of the genomic range.

        Might be useful for saving to json file format.

        Returns:
            dict: Dictionary representation of the genomic range.
        """
        return {'chrom': self.chrom,
                'start': self.start,
                'end': self.end,
                'strand': self.strand,
                'genome': self.genome,
                'sequence_file_path': self.sequence_file_path.to_dict(),
                'name': self.name,
                'relations': {key: value.to_dict()
                              for key, value in self.relations.items()},
                'relations_class': (next(iter(self.relations.values())).__class__.__name__
                                    if len(self.relations) > 0 else None)}

    @classmethod
    def from_dict(cls, dict_):
        """Gets genomic range from its dictionary representation.

        Args:
            dict_ (dict): A dictionary generated with
                `GenomicRange.to_dict`.

        Returns:
            GenomicRange: Restored genomic range.

        Raises:
            ValueError: An error in case key entry *relation_class* is not one
                of `BaseGenomicRangesList`, `GenomicRangesList` or None.
            ValueError: An error in case *dict_* does not contain any
                one of necessary keys.
        """
        relations_class_name = dict_.get('relations_class', 'GenomicRangesList')
        if relations_class_name == 'BaseGenomicRangesList':
            relations_class = BaseGenomicRangesList
        elif relations_class_name == 'GenomicRangesList':
            relations_class = GenomicRangesList
        elif relations_class_name is None:
            relations_class = GenomicRangesList
        else:
            raise ValueError("Provided relations_class value is not one of "
                             "'BaseGenomicRangesList', 'GenomicRangesList', None.")
        parser_dict = {'chrom': lambda i: i.get('chrom'),
                       'start': lambda i: i.get('start'),
                       'end': lambda i: i.get('end'),
                       'strand': lambda i: i.get('strand'),
                       'genome': lambda i: i.get('genome'),
                       'name': lambda i: i.get('name'),
                       'sequence_file_path': lambda i: SequencePath.from_dict(i.get('sequence_file_path', {'path': None})),
                       'relations': lambda i: {key: relations_class.from_dict(value)
                                               for key, value in i.get('relations', {}).items()}}
        for key in parser_dict:
            if key not in dict_:
                raise ValueError(f"Provided dict_ does not contain {key} key.")
        return cls(**{key: func(dict_) for key, func in parser_dict.items()})

    @property
    def sequence_file_path(self):
        return self._sequence_file_path

    @sequence_file_path.setter
    def sequence_file_path(self, value):
        self._sequence_file_path = SequencePath(value)

    def distance(self, other):
        """Calculates distance between two genomic ranges.

        In case chromosomes and genomes are consistent, returns distance
        between two genomic ranges in strand-nonspecific way.

        Args:
            other (BaseGenomicRange or children): The other genomic range
                to calculate distance to *self*.

        Raises:
            InconsistentGenomesError: An error in case of inconsistent genomes.
            InconsistentChromosomesError: An error in case of inconsistent chromosomes.
        """
        if hasattr(other, 'genome') and self.genome != other.genome:
            raise InconsistentGenomesError(self.genome, other.genome)
        return super().distance(other)

    def align_blast(self, other, program='blastn', task='blastn', **kwargs):
        """Aligns two genomic ranges with BLAST.

        If both genomic ranges have their sequences
        stored in separate files, aligns them with
        standalone blastn and specified parameters.
        Alignment file format is "7 std score".
        Standalone blast package should be installed
        manually.

        Args:
            other (GenomicRange): Genomic range to
                align self with.
            program (str, default "blastn"): A program name
                from Standalone BLAST to use.
            task (str, default "blastn"): Which task from a BLAST
                program to use.
            kwargs (dict): Any command line arguments
                to blastn in the following scheme:
                ``{'flag': 'value'}``.

        Returns:
            GenomicRangesAlignment: Aligned range pair, which associates
            self and other genomic ranges and their alignment.

        Raises:
            FilePathNotSpecifiedError: In case one or both
                genomic ranges do not have separate
                sequence files.
        """
        for grange in (self, other):
            grange.sequence_file_path.check_correct()

        if program not in ['blastn', 'blastp']:
            raise ValueError("`program` attribute is not one of "
                             "['blastn', 'blastp']")

        command = [program,
                   '-task', task,
                   '-query', self.sequence_file_path,
                   '-subject', other.sequence_file_path,
                   '-outfmt', "7 std score"]
        add_args = sum((['-' + str(key), str(value)]
                        for key, value in kwargs.items()),
                       [])
        with Popen(command + add_args, stdout=PIPE, stderr=PIPE) as proc:
            stream = TextIOWrapper(proc.stdout)
            alignment = GenomicRangesAlignment.from_file_blast(stream,
                                                               self,
                                                               other)
            stream.close()
        return alignment


class GenomicRangesAlignment(Alignment):
    __slots__ = ('qrange', 'srange', 'relative')

    def __init__(self,
                 HSPs,
                 qrange,
                 srange,
                 qlen=None,
                 slen=None,
                 filtered_HSPs=None,
                 filtered=None,
                 relative=True):
        """Inits Alignment class.

        Args:
            HSPs (list): list of HSP instances;
            qrange (GenomicRange): query genomic range;
            srange (GenomicRange): subject genomic range;
            qlen (int): query sequence length
                (default: None);
            slen (int): subject sequence length
                (default: None);
            filtered_HSPs (list): list of HSP instances
                that were filtered by some criterion
               (default: None).

        """
        super().__init__(HSPs, qlen, slen, filtered_HSPs, filtered)
        self.qrange = qrange
        self.srange = srange
        self.relative = relative

    def __str__(self):
        return super().__str__()

    def __repr__(self):
        return super().__str__()

    def __eq__(self, other):
        if len(self._all_HSPs) != len(other._all_HSPs):
            return False
        if self.qlen != other.qlen or self.slen != other.slen:
            return False
        if self.qrange != other.qrange or self.srange != other.srange:
            return False
        if self.relative != other.relative:
            return False
        return all([i == k for i, k in zip(self._all_HSPs, other._all_HSPs)])

    @property
    def chain_class(self):
        return GenomicRangesAlignmentChain

    def to_dict(self):
        dict_ = super().to_dict()
        dict_.update({'qrange': self.qrange.to_dict(),
                      'srange': self.srange.to_dict(),
                      'relative': self.relative})
        return dict_

    @classmethod
    def from_dict(cls, dict_):
        keys = ['HSPs',
                'qlen',
                'slen',
                'filtered',
                'filtered_HSPs',
                'qrange',
                'srange',
                'relative']
        functions = [lambda i: [HSPVertex.from_dict(item)
                                for item in i.get('HSPs')],
                     lambda i: i.get('qlen'),
                     lambda i: i.get('slen'),
                     lambda i: i.get('filtered'),
                     lambda i: None if i.get('filtered_HSPs') is None else [HSPVertex.from_dict(item)
                                                                            for item in i.get('HSPs')],
                     lambda i: GenomicRange.from_dict(i.get('qrange')),
                     lambda i: GenomicRange.from_dict(i.get('srange')),
                     lambda i: i.get('relative')]
        kwargs = {key: function(dict_)
                  for key, function in zip(keys, functions)}
        return cls(**kwargs)

    def to_genomic(self):
        if not self.relative:
            return
        # query sequence
        if self.qrange.strand != '-':
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.qstart = hsp.qstart + self.qrange.start
                hsp.qend = hsp.qend + self.qrange.start
        else:
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.qstart = self.qrange.end - hsp.qstart
                hsp.qend = self.qrange.end - hsp.qend
        # subject sequence
        if self.srange.strand != '-':
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.sstart = hsp.sstart + self.srange.start
                hsp.send = hsp.send + self.srange.start
        else:
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.sstart = self.srange.end - hsp.sstart
                hsp.send = self.srange.end - hsp.send
        self.relative = False

    def to_relative(self):
        if self.relative:
            return
        if self.qrange.strand != '-':
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.qstart = hsp.qstart - self.qrange.start
                hsp.qend = hsp.qend - self.qrange.start
        else:
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.qstart = self.qrange.end - hsp.qstart
                hsp.qend = self.qrange.end - hsp.qend
        # subject sequence
        if self.srange.strand != '-':
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.sstart = hsp.sstart - self.srange.start
                hsp.send = hsp.send - self.srange.start
        else:
            for hsp in chain(self._all_HSPs, self.HSPs):
                hsp.sstart = self.srange.end - hsp.sstart
                hsp.send = self.srange.end - hsp.send
        self.relative = True

    def cut_coordinates(self, qleft=None, qright=None, sleft=None, sright=None):
        hsps = self._cut_hsps(qleft=qleft, qright=qright, sleft=sleft, sright=sright)
        return GenomicRangesAlignment(hsps,
                                      qrange=self.qrange,
                                      srange=self.srange,
                                      relative=self.relative)

    @classmethod
    def from_file_blast(cls,
                        file_object,
                        qrange,
                        srange,
                        start_type=1,
                        end_type='inclusive'):
        return super().from_file_blast(file_object,
                                       start_type=start_type,
                                       end_type=end_type,
                                       qrange=qrange,
                                       srange=srange)


class GenomicRangesAlignmentChain(AlignmentChain):

    def __init__(self, HSPs, alignment, score=None):
        super().__init__(HSPs, alignment, score)

    @classmethod
    def from_dict(cls, dict_):
        keys = ['HSPs',
                'alignment',
                'score']
        functions = [lambda i: [HSPVertex.from_dict(item)
                                for item in i.get('HSPs', [])],
                     lambda i: GenomicRangesAlignment.from_dict(i.get('alignment')),
                     lambda i: i.get('score')]
        kwargs = {key: function(dict_)
                  for key, function in zip(keys, functions)}
        return cls(**kwargs)

    def _prepare_side(self, side='q'):
        if side == 'q':
            chrom = self.alignment.qrange.chrom
            chromStart = min(min(self.HSPs, key=lambda i: i.qstart).qstart,
                             min(self.HSPs, key=lambda i: i.qend).qend)
            chromEnd = max(max(self.HSPs, key=lambda i: i.qstart).qstart,
                           max(self.HSPs, key=lambda i: i.qend).qend)
            name = self.alignment.qrange.name
            strand = nxor_strands(self.HSPs[0].qstrand,
                                  self.alignment.qrange.strand)
            blockSizes = [abs(hsp.qend - hsp.qstart) for hsp in self.HSPs]
            if strand == "+":
                blockStarts = [hsp.qstart - chromStart for hsp in self.HSPs]
            else:
                blockStarts = [hsp.qend - chromStart for hsp in self.HSPs]
        elif side == 's':
            chrom = self.alignment.srange.chrom
            chromStart = min(min(self.HSPs, key=lambda i: i.sstart).sstart,
                             min(self.HSPs, key=lambda i: i.send).send)
            chromEnd = max(max(self.HSPs, key=lambda i: i.sstart).sstart,
                           max(self.HSPs, key=lambda i: i.send).send)
            name = self.alignment.srange.name
            strand = nxor_strands(self.HSPs[0].sstrand,
                                  self.alignment.srange.strand)
            blockSizes = [abs(hsp.send - hsp.sstart) for hsp in self.HSPs]
            if strand == "+":
                blockStarts = [hsp.sstart - chromStart for hsp in self.HSPs]
            else:
                blockStarts = [hsp.send - chromStart for hsp in self.HSPs]
        else:
            raise ValueError(f'side argument {side} is not one of ["q", "s"]')
        score = self.score
        thickStart = chromStart
        thickEnd = chromEnd
        itemRgb = 0
        blockCount = len(self.HSPs)
        side = [chrom,
                chromStart,
                chromEnd,
                name,
                score,
                strand,
                thickStart,
                thickEnd,
                itemRgb,
                blockCount,
                blockSizes,
                blockStarts]
        return side

    def to_bed12(self, mode='list'):
        """
        Turns alignemnt chain into BED12 representation of both
        query side and subject side.

        Args:
            mode (str): how to return representation. If 'list' then
            list of lists, if 'str' then complete BED12 record.

        Returns:
            (list of two lists or two strs): BED12 representation.
        """
        if len(self.HSPs) == 0:
            raise ValueError('No HSPs were found in the AlignmentChain.')
        if self.alignment.relative:
            raise ValueError('Alignment coordinates are not translated to genomic coordinates. '\
                             'Use GenomicRangesAlignmentChain.alignment.to_genomic().')
        q_side = self._prepare_side(side='q')
        s_side = self._prepare_side(side='s')
        if mode == 'list':
            return [q_side, s_side]
        elif mode == 'str':
            return ["\t".join([','.join(str(value) for value in item) if isinstance(item, list) else str(item)
                               for item in q_side]),
                    "\t".join([','.join(str(value) for value in item) if isinstance(item, list) else str(item)
                               for item in s_side])]
        else:
            raise ValueError('mode not one of ["list", "str"].')


ChromosomeLocation = namedtuple('ChromosomeLocation', 'size start')


class FastaSeqFile:
    """Represents genomic fasta file.

    Should locate to the fasta file
    where the genome is placed. Chromosome
    names must be placed in fasta record
    header and separated from other header
    information with whitespace character.

    Attributes:
        sequence_file_path (Path): path to the file.
        chromsizes (dict): dictionary of
            ChromosomseLocation instances.
        _file_obj: opened file under the
            sequence_file_path.
        line_length: line length of fasta sequence.
    """

    def __init__(self, sequence_file_path):
        """Initializes FastaSeqFile instance.

        Args:
            sequence_file_path (str, Path): path to the file.
        """
        self._sequence_file_path = SequencePath(sequence_file_path)
        self._chromsizes = dict()
        self._file_obj = None
        self._line_length = None

    def __repr__(self):
        """Returns representation of FastaSeqFile object."""
        return f"FastaSeqFile('{self.sequence_file_path}')"

    def __str__(self):
        """Returns string representation."""
        return f"FastaSeqFile('{self.sequence_file_path}')"

    def __eq__(self, other):
        """Equality magic method."""
        return self.sequence_file_path == other.sequence_file_path

    def to_dict(self):
        """Returns dict representation."""
        return {'sequence_file_path': self.sequence_file_path.to_dict()}

    @classmethod
    def from_dict(cls, dict_):
        """Restores FastaSeqFile instance from its dict representation."""
        keys = ['sequence_file_path']
        try:
            kwargs = {key: SequencePath.from_dict(dict_[key]) for key in keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of nececcary keys {keys}.")
        return cls(**kwargs)

    @property
    def sequence_file_path(self):
        """`sequence_file_path` getter"""
        return self._sequence_file_path

    @sequence_file_path.setter
    def sequence_file_path(self, value):
        """`sequence_file_path` setter"""
        self._sequence_file_path = SequencePath(value)

    @property
    def chromsizes(self):
        """Property representing genomic chromosomes.

        Calculates chromosome names, their sizes and
        sequence start positions in the genome file.

        Returns:
            (dict) dictionary of ChromosomeLocation
                instances in the following scheme:
                {chromosome_name: ChromosomeLocation(size, start)}

        Raises:
            FilePathNotSpecifiedError in case the location of
                sequence file is not specified.
        """
        if not self._chromsizes:

            self.sequence_file_path.check_correct()

            with open(self.sequence_file_path, 'r') as infile:
                chrom, size, start = None, None, None
                line = infile.readline()
                while line:
                    if line.startswith(">"):
                        self._chromsizes[chrom] = ChromosomeLocation(size,
                                                                     start)
                        chrom = line.lstrip(">").rstrip().split()[0]
                        size = 0
                        start = None
                    else:
                        if start is None:
                            start = infile.tell() - len(line)
                        size += len(line.strip())
                    line = infile.readline()
                self._chromsizes[chrom] = ChromosomeLocation(size, start)
                del self._chromsizes[None]
        return self._chromsizes

    @property
    def line_length(self):
        """`line_length` getter.

        Returns:
            (int) length of sequence line in fasta file.
        """
        if self._line_length is None:
            with open(self.sequence_file_path, 'r') as infile:
                line = infile.readline()
                while line.startswith(">"):
                    line = infile.readline()
                self._line_length = len(line.strip())
        return self._line_length

    def __enter__(self):
        """Context manager for opening sequence file."""
        self.sequence_file_path.check_correct()
        self._file_obj = open(self.sequence_file_path, 'r')

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager for closing sequence file."""
        self._file_obj.close()
        self._file_obj = None

    def locate_coordinate(self, chrom_start, coord):
        """Puts file pointer to the beginning of the genomic range.

        Args:
            chrom_start (int): byte number of chromosome start.
            coord (int): coordinate in the corresponding chromosome.

        Returns:
            None
        """
        self._file_obj.seek(chrom_start + coord + coord // self.line_length)

    def get_fasta_by_coord(self, grange, outfile, chromsizes=None, linewidth=60):
        """Writes sequence of the genomic range in fasta file.

        Writes full fasta record of the given genomic range
        to the file. The header is made in the following
        scheme:
        ">" + grange.sequence_header + "\n"

        Args:
            grange (GenomicRange): genomic range which fasta
                sequence write to file.
            outfile (fileobj): file object to which write fasta
                sequence.
            chromsizes (dict): a dict of ChromosomeLocation tuples
                derived from self.chromsizes. If None,
                estimates itself, may affect performance (default: None).
            linewidth (int): line length of fasta sequence to write
                (default: 60).

        Returns:
            None.

        Raises:
            ChromosomeNotFoundError if the chromosome of
                the given genomic range is not presented
                in FastaSeqFile instance.
        """
        if chromsizes is None:
            chromsizes = self.chromsizes
        try:
            chrom_size, chrom_start = chromsizes[grange.chrom]
        except KeyError:
            raise ChromosomeNotFoundError(grange.chrom,
                                          self.sequence_file_path,
                                          chromsizes.keys())
        if grange.end > chrom_size:
            raise GenomicCoordinatesError(grange,
                                          chromsizes[grange.chrom])
        self.locate_coordinate(chrom_start, grange.start)
        if grange.strand == '-':
            reverse = True
            complement = True
        else:
            reverse = False
            complement = False
        outfile.write(">" + grange.sequence_header + "\n")
        lines_before_start = grange.start // self.line_length
        lines_before_end = grange.end // self.line_length
        lines_inside = lines_before_end - lines_before_start
        total = grange.end - grange.start + lines_inside
        fasta_getter(self._file_obj,
                     outfile,
                     total=total,
                     linewidth=linewidth,
                     reverse=reverse,
                     complement=complement)


# A constant defining how many genomic ranges can be printed
# inside the genomic ranges list via repr of str.
MAX_ROWS = 10


class BaseGenomicRangesList(SortedKeyList):
    """Basic sorted mutable collection of genomic ranges.

    Subclass of `SortedKeyList`. Provides sorting of
    `BaseGenomicRange` and its children. Requiers any item
    to have 3 attributes for sorting: `chrom`, `start`, `end`.

    Attributes:

    Basic methods:
        __init__
        __repr__
        __str__
        __eq__

    Dictionary interoperability:
        to_dict
        from_dict

    Genomic methods:
        merge
        close_merge
        inter_ranges
        basic_annotation_parser

    Properties:
        init_args

    Static methods:
        _genomic_key_func

    Class methods:
        from_dict
        basic_annotation_parser
    """

    @staticmethod
    def _genomic_key_func(item):
        """Returns item attributes as key for sorting.

        Returns:
            tuple: chrom, start, end.
        """
        return item.chrom, item.start, item.end

    def __new__(cls, collection=[], *args, **kwargs):
        """Makes new instance of BaseGenomicRangesList."""
        instance = super().__new__(cls,
                                   iterable=collection,
                                   key=cls._genomic_key_func)
        return instance

    def __init__(self, collection):
        """Initializes BaseGenomicRangesList instance.

        Args:
            collection (iterable): any collection containing
                BaseGenomicRange instances.

        Returns:
            None
        """
        super().__init__(iterable=collection, key=self._genomic_key_func)

    @property
    def init_args(self):
        """Returns tuple of __init__ arguments."""
        return ('collection',)

    def __repr__(self):
        """Returns code representation of the list of genomic ranges.

        A number of genomic ranges listed is controlled by the constant
        `ortho2align.genomicranges.MAX_ROWS`. Change it to `None` to dismiss
        limitation.
        """
        max_rows = len(self) if MAX_ROWS is None else MAX_ROWS
        granges = ", ".join([repr(grange) for grange in self][:max_rows])
        appendix = ", ..." if max_rows < len(self) else ""
        return f"BaseGenomicRangesList([{granges}{appendix}])"

    def __str__(self):
        """Returns tab-separated table of containing genomic ranges.

        A number of genomic ranges listed is controlled by the constant
        `ortho2align.genomicranges.MAX_ROWS`. Change it to `None` to dismiss
        limitation.
        """
        max_rows = len(self) if MAX_ROWS is None else MAX_ROWS
        granges = "\n".join([str(grange) for grange in self][:max_rows])
        appendix = "\n..." if max_rows < len(self) else ""
        return f"BaseGenomicRangesList of {len(self)} genomic ranges:\n"  \
               f"{granges}{appendix}\n" \
               f"{min(max_rows, len(self))} genomic ranges are shown."

    def __eq__(self, other):
        """Checks for equality of two lists of genomic ranges.

        First checks if lengths of two lists are equal. If `True`,
        checks equality of all genomic ranges in sorted order.
        """
        if len(self) != len(other):
            return False
        return all([self_grange == other_grange
                    for self_grange, other_grange in zip(self, other)])

    def to_dict(self):
        """Returns dictionary representation of the list of genomic ranges.

        Returns:
            dict: 'collection' key stands for the list of dictionarized
            genomic ranges, whereas 'grange_class' key stands for the name
            of the class of the first genomic range or None in case of empty list.
        """
        return {'collection': [grange.to_dict() for grange in self],
                'grange_class': self[0].__class__.__name__ if len(self) > 0 else None}

    @classmethod
    def from_dict(cls, dict_):
        """Restores list of genomic ranges from its dictionary representation.

        Args:
            dict_ (dict): a dictionary generated with `to_dict` method. Must
                contain 2 keys: 'collection' for the list of the dictionarized
                genomic ranges and 'grange_class' for the name of the class
                of genomic ranges. 'grange_class' must be one of two:
                `BaseGenomicRange` or `GenomicRange`.

        Returns:
            BaseGenomicRangesList: restored list of genomic ranges.

        Raises:
            ValueError: an error in case there is no one of the essential
            keywords in the `dict_` or 'grange_class' is not one of `BaseGenomicRange`,
            `GenomicRange`.
        """
        grange_classname = dict_.get('grange_class', 'GenomicRange')
        if grange_classname == 'BaseGenomicRange':
            grange_class = BaseGenomicRange
        elif grange_classname == 'GenomicRange':
            grange_class = GenomicRange
        else:
            raise ValueError(f"Provided `grange_class` value {grange_classname} "
                             f"is not one of acceptable: [`BaseGenomicRange`, `GenomicRange`].")
        parser_dict = {'collection': lambda i: [grange_class.from_dict(item)
                                                for item in i.get('collection')]}
        for key in parser_dict:
            if key not in dict_:
                raise ValueError('There is no key `{key}` in the provided dictionary.')
        return cls(**{key: func(dict_) for key, func in parser_dict.items()})

    def merge(self, distance=0, verbose=False):
        """Merges all genomic ranges in the list on the given distance.

        In case two genomic ranges are of the same chromosome and
        are not further from each other than the provided distance value,
        they are merged together in strand-nonspecific way. This method
        will be fast in case of sparse list of genomic ranges (i.e. when
        there is no many close genomic ranges).

        Args:
            distance (int): max distance between two genomic ranges
                that can be merged together (default: 0).
            verbose (bool): whether to show progress via `tqdm`
                progress bar or not (default: False).

        Returns:
            BaseGenomicRangesList: new list with merged genomic ranges.
        """
        new_list = list()
        if len(self) != 0:
            new_list.append(self[0])
            for grange in tqdm(self[1:],
                               desc='Merging genomic ranges',
                               unit='grange',
                               disable=not verbose):
                try:
                    if abs(new_list[-1].distance(grange)) <= distance:
                        old_range = new_list.pop()
                        new_list.append(old_range.merge(grange))
                    else:
                        new_list.append(grange)
                except GenomicException:
                    new_list.append(grange)

        used_args = {'collection', }
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(new_list, **kwargs)

    def close_merge(self, distance=0, verbose=False):
        """Merges all genomic ranges in the list on the given distance.

        In case two genomic ranges are of the same chromosome and
        are not further from each other than the provided distance value,
        they are merged together in strand-nonspecific way. This method
        will be fast in case of dense list of genomic ranges (i.e. when
        there are many close genomic ranges).

        Args:
            distance (int): max distance between two genomic ranges
                that can be merged together (default: 0).
            verbose (bool): whether to show progress via `tqdm`
                progress bar or not (default: False).

        Returns:
            BaseGenomicRangesList: new list with merged genomic ranges.
        """
        new_list = list()
        if len(self) != 0:
            new_list.append(self[0])
            prev_grange = self[0]
            for grange in tqdm(self[1:],
                               desc='Merging genomic ranges',
                               unit='grange',
                               disable=not verbose):
                try:
                    if abs(prev_grange.distance(grange)) > distance:
                        old_range = new_list.pop()
                        new_list.append(old_range.merge(prev_grange))
                        new_list.append(grange)
                except GenomicException:
                    old_range = new_list.pop()
                    new_list.append(old_range.merge(prev_grange))
                    new_list.append(grange)
                prev_grange = grange
            # in case there were no last merge in loop
            old_range = new_list.pop()
            new_list.append(old_range.merge(prev_grange))

        used_args = {'collection', }
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(new_list, **kwargs)

    def inter_ranges(self, distance=0, verbose=False):
        """Returns ranges between merged ranges.

        Inter-ranges are genomic ranges between presented
        genomic ranges (i.e. they start at end coordinates
        and end at start coordinates of merged genomic ranges).
        First, it merges all genomic ranges at the given distance,
        next inverts genomic ranges. Strand of inter-ranges is `"."`.

        Args:
            distance (int): max distance between two genomic ranges
                that can be merged together (default: 0).
            verbose (bool): whether to show progress via `tqdm`
                progress bar or not (default: False).

        Returns:
            BaseGenomicRangesList: a new list with inter-ranges.
        """
        merged = self.merge(distance)
        inverted_granges = list()
        if len(merged) != 0:
            used_grange_args = {'chrom', 'start', 'end', 'strand'}
            for i in tqdm(range(len(merged) - 1),
                          desc='Finding intergenomic ranges',
                          unit='grange',
                          disable=not verbose):
                if merged[i].chrom == merged[i + 1].chrom:
                    new_chrom = merged[i].chrom
                    new_start = merged[i].end
                    new_end = merged[i + 1].start
                    new_strand = "."
                    grange_kwargs = {attr: getattr(merged[i], attr)
                                     for attr in set(merged[i].init_args) - used_grange_args}
                    inverted_granges.append(merged[i].__class__(new_chrom,
                                                                new_start,
                                                                new_end,
                                                                new_strand,
                                                                **grange_kwargs))

        used_args = {'collection', }
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(inverted_granges, **kwargs)

    def to_bed6(self, fileobj):
        """Writes contents in bed6 format to a file object.

        Args:
            fileobj (file-like): file-like object for output.
        """
        for grange in self:
            fileobj.write('\t'.join([str(grange.chrom),
                                     str(grange.start),
                                     str(grange.end),
                                     str(grange.name),
                                     '.',
                                     str(grange.strand)]) + '\n')

    @classmethod
    def basic_annotation_parser(cls, fileobj, **kwargs):
        raise NotImplementedError


class GenomicRangesList(BaseGenomicRangesList):
    """Represents a list of `GenomicRange` instances.

    Attributes:
        sequence_file_path (SequencePath): A path to the corresponding
            genome file.
        _sequence_file_path (SequencePath): A placeholder for the
            *sequence_file_path* property.
        sequence_file (FastaSeqFile): `FastaSeqFile` instance with the
            corresponding genome file.
        name_mapping (dict of lists): A dictionary with keys
            as genomic ranges names and values as lists of
            corresponding GenomicRange instances.
        init_args (tuple): Essential *__init__* arguments.

    Class attributes:
        column_names (dict of lists): Column names of selected
            file formats of genomic annotations.
        name_columns (dict of strs): Names of columns containing
            information about genomic range names for selected
            file formats of genomic annotations.
        comments (dict of strs): Comment characters for selected file
            formats of genomic annotations.
        name_patterns (dict of r-strs): Regexp patterns to extract
            genomic range from the column containing the name for
            selected file formats of genomic annotations.
        seps (dict of strs): Field separator characters for selected
            file formats of genomic annotations.
        fileformats (list): List of selected file formats of genomic
            annotations (gtf, gff, bed3, bed6, bed12, None &mdash; user-specified).
        dtypes (dict of lists): Dict of lists of callables representing
            dtypes of columns for selected file formats of genomic
            annotations.
        start_types (dict of ints): Maps file format to genomic range
            coordinate start type (0 or 1).
        end_types (dict of strs): Maps file format to genomic range
            coordinate end type ('inclusive' or 'exclusive').
    """
    column_names = {'gff': ['chrom',
                            'source',
                            'type',
                            'start',
                            'end',
                            'smth1',
                            'strand',
                            'smth2',
                            'data'],
                    'gtf': ['chrom',
                            'source',
                            'type',
                            'start',
                            'end',
                            'smth1',
                            'strand',
                            'smth2',
                            'data'],
                    'bed3': ['chrom',
                             'start',
                             'end'],
                    'bed6': ['chrom',
                             'start',
                             'end',
                             'data',
                             'score',
                             'strand'],
                    'bed12': ['chrom',
                              'start',
                              'end',
                              'data',
                              'score',
                              'strand',
                              'thickstart',
                              'thickend',
                              'itemrgb',
                              'blockcount',
                              'blocksizes',
                              'blockstarts']}
    name_columns = {'gtf': 'data',
                    'gff': 'data',
                    'bed3': None,
                    'bed6': 'data',
                    'bed12': 'data'}
    comments = {'gtf': '#',
                'gff': '#',
                'bed3': '#',
                'bed6': '#',
                'bed12': '#'}
    name_patterns = {'gff': r'GeneID:(\d+)',
                     'gtf': r'gene_id\s\"(\w*)\"',
                     'bed3': r'(.*)',
                     'bed6': r'(.*)',
                     'bed12': r'(.*)'}
    seps = {'gff': '\t',
            'gtf': '\t',
            'bed3': '\t',
            'bed6': '\t',
            'bed12': '\t'}
    fileformats = ['gtf',
                   'gff',
                   'bed3',
                   'bed6',
                   'bed12',
                   None]
    dtypes = {'gff': [str,
                      str,
                      str,
                      int,
                      int,
                      str,
                      str,
                      str,
                      str],
              'gtf': [str,
                      str,
                      str,
                      int,
                      int,
                      str,
                      str,
                      str,
                      str],
              'bed3': [str,
                       int,
                       int],
              'bed6': [str,
                       int,
                       int,
                       str,
                       numberize,
                       str],
              'bed12': [str,
                        int,
                        int,
                        str,
                        numberize,
                        str,
                        int,
                        int,
                        str,
                        int,
                        lambda x: x.split(","),
                        lambda x: x.split(",")]}
    start_types = {'gtf': 1,
                   'gff': 1,
                   'bed3': 0,
                   'bed6': 0,
                   'bed12': 0}
    end_types = {'gtf': 'inclusive',
                 'gff': 'inclusive',
                 'bed3': 'exclusive',
                 'bed6': 'exclusive',
                 'bed12': 'exclusive'}

    def __init__(self, collection, sequence_file_path):
        """Initializes GenomicRangesList instance.

        Args:
            collection (iterable of `GenomicRange`): A collection of 
                `GenomicRange` instances.
            sequence_file_path (str, SequencePath): A path to the
                corresponding genome file in fasta format.

        Returns:
            None
        """
        super().__init__(collection=collection)
        self._sequence_file_path = (sequence_file_path
                                    if isinstance(sequence_file_path,
                                                  SequencePath)
                                    else SequencePath(sequence_file_path))
        self.sequence_file = FastaSeqFile(self.sequence_file_path)
        self._name_mapping = defaultdict(list)

    @property
    def init_args(self):
        """Returns a tuple of essential *__init__* arguments."""
        return ('collection', 'sequence_file_path')

    def __repr__(self):
        """Returns code representation of the list of genomic ranges.

        A number of genomic ranges listed is controlled by the constant
        `ortho2align.genomicranges.MAX_ROWS`. Change it to `None` to dismiss
        limitation.
        """
        max_rows = len(self) if MAX_ROWS is None else MAX_ROWS
        granges = ", ".join([repr(grange) for grange in self][:max_rows])
        appendix = ", ..." if max_rows < len(self) else ""
        return f"GenomicRangesList({granges}{appendix}, {self.sequence_file_path})"

    def __str__(self):
        """Returns tab-separated table of containing genomic ranges.

        A number of genomic ranges listed is controlled by the constant
        `ortho2align.genomicranges.MAX_ROWS`. Change it to `None` to dismiss
        limitation.
        """
        max_rows = len(self) if MAX_ROWS is None else MAX_ROWS
        granges = "\n".join([str(grange) for grange in self][:max_rows])
        appendix = "\n..." if max_rows < len(self) else ""
        return f"GenomicRangesList of {len(self)} genomic ranges from {self.sequence_file_path}:\n"  \
               f"{granges}{appendix}\n" \
               f"{min(max_rows, len(self))} genomic ranges are shown."

    def __eq__(self, other):
        """Checks for equality of two lists of genomic ranges."""
        if len(self) != len(other):
            return False
        if self.sequence_file_path != other.sequence_file_path:
            return False
        return all([self_grange == other_grange
                    for self_grange, other_grange in zip(self, other)])

    def to_dict(self):
        """Returns dictionary representation of the list of genomic ranges.

        Returns:
            dict: 'collection' key stands for the list of dictionarized
            genomic ranges, 'sequence_file_path' key stands for dictionary
            representation of corresponding attribute, grange_class' key stands
            for the name of the class of the first genomic range or `None` in case
            of empty list.
        """
        return {'collection': [grange.to_dict() for grange in self],
                'sequence_file_path': self.sequence_file_path.to_dict(),
                'grange_class': self[0].__class__.__name__ if len(self) > 0 else None}

    @classmethod
    def from_dict(cls, dict_):
        """Restores list of genomic ranges from its dictionary representation.

        Args:
            dict_ (dict): A dictionary generated with `to_dict` method. Must
                contain 3 keys: 'collection' for the list of the dictionarized
                genomic ranges, 'grange_class' for the name of the class
                of genomic ranges and 'sequence_file_path' for the path
                to the corresponding genome file in fasta format.
                'grange_class' must be one of two: `BaseGenomicRange` or `GenomicRange`.

        Returns:
            GenomicRangesList: Restored list of genomic ranges.

        Raises:
            ValueError: An error in case there is no one of the essential
                keywords in the `dict_` or 'grange_class' is not one of
                `BaseGenomicRange`, `GenomicRange`.
        """
        grange_classname = dict_.get('grange_class', 'GenomicRange')
        if grange_classname == 'BaseGenomicRange':
            grange_class = BaseGenomicRange
        elif grange_classname == 'GenomicRange':
            grange_class = GenomicRange
        else:
            raise ValueError(f"Provided `grange_class` value {grange_classname} "
                             f"is not one of acceptable: [`BaseGenomicRange`, `GenomicRange`].")
        parser_dict = {'collection': lambda i: [grange_class.from_dict(item)
                                                for item in i.get('collection')],
                       'sequence_file_path': lambda i: SequencePath.from_dict(i.get('sequence_file_path',
                                                                                    {'path': None}))}
        for key in parser_dict:
            if key not in dict_:
                raise ValueError('There is no key `{key}` in the provided dictionary.')
        return cls(**{key: func(dict_) for key, func in parser_dict.items()})

    @property
    def sequence_file_path(self):
        return self._sequence_file_path

    @sequence_file_path.setter
    def sequence_file_path(self, value):
        self._sequence_file_path = (value
                                    if isinstance(value,
                                                  SequencePath)
                                    else SequencePath(value))

    @property
    def name_mapping(self):
        if not self._name_mapping:
            for grange in self:
                self._name_mapping[grange.name].append(grange)
        return self._name_mapping

    def reset_name_mapping(self):
        """Resets name mapping.

        Must be called after changes in collection
        and before calling `GenomicRangesList.name_mapping` property.
        """
        self._name_mapping = defaultdict(list)

    def get_neighbours(self, other, distance=0, relation='neighbours', verbose=False):
        """Finds all neighbours of genomic ranges in *self* and assigns them to *relations*.

        For every genomic range in *self* finds neighbours in *other*
        at given *distance* and assigns a list of them to the specified
        *relation*. Can be used in case genomic ranges in *self* are instances
        of `GenomicRange`.

        Args:
            other (BaseGenomicRangesList and children):
                The other genomic ranges list to find neighbours in.
            distance (int, default 0): Distance to find neighbours at.
            relation (str, default 'neighbours'): Name of the relation
                to assign neighbours in.
            verbose (bool, default False): Whether to report progress
                via tqdm progress bar or not.
        """
        for grange in tqdm(self,
                           desc='Getting neighbours',
                           unit='grange',
                           disable=not verbose):
            grange.relations[relation] = grange.find_neighbours(other,
                                                                distance=distance)

    def flank(self, distance=0, check_boundaries=True, chromsizes=None, verbose=False):
        """Flanks genomic ranges concerning chromosome boundaries.

        After flanking (see `BaseGenomicRange.flank`) optionally
        checks whether flanked genomic ranges are out of chromosome
        boundaries. Start coordinate sets to zero in case it's negative.
        End coordinate sets to chromosome size in case it's bigger than
        latter. Requires *chromsizes* dict to be passed or
        `GenomicRangesList.sequence_file_path` to be a real path
        to the corresponding genome file. Indexing genome in the latter
        case might take a while.

        Args:
            distance (int, default 0): A distance to flank genomic ranges.
            check_boundaries (bool, default True): Whether to check
                for violation of chromosome boundaries by flanked genomic
                ranges or not.
            chromsizes (dict of ChromosomeLocation, default None): A dictionary
                with chromosome names as keys and `ChromosomeLocation` named tuples
                as values. If *None*, calculates itself from the *sequence_file* via
                `FastaSeqFile.chromsizes`.
            verbose (bool, default False): Whether to report progress via
                tqdm progress ar or not.

        Returns:
            self.__class__: A new genomic ranges list with flanked genomic ranges.
        """
        new_list = list()
        if chromsizes is None:
            chromsizes = self.sequence_file.chromsizes
        for grange in tqdm(self,
                           desc='Flanking genomic ranges',
                           unit='grange',
                           disable=not verbose):
            new_grange = grange.flank(distance)
            if check_boundaries:
                new_grange.start = max(0, new_grange.start)
                new_grange.end = min(chromsizes[new_grange.chrom].size,
                                     new_grange.end)
            new_list.append(new_grange)
        used_args = {'collection', }
        kwargs = {attr: getattr(self, attr)
                  for attr in set(self.init_args) - used_args}
        return self.__class__(new_list, **kwargs)

    def get_fasta(self, outfileprefix, outdir=SequencePath('.'), mode='split', chromsizes=None, verbose=False):
        """Extracts sequences of genomic ranges.

        Can extract sequences of all genomic ranges in *'bulk'* mode
        into single file named ``*outfileprefix* + '.fasta'``
        or each sequence into separate files named
        ``*outfileprefix* + genomic_range.name + '.fasta'`` in *'split'* mode .

        Args:
            outfileprefix (str, Path, SequencePath): Prefix of the
                output file(s).
            outdir (str, Path, SequencePath, default SequencePath('.')): Output directory.
            mode (str, default 'split'): The mode to extract sequences,
                one of 'split', 'bulk'.
            chromsizes (dict, default None): A dict of `ChromosomeLocation`
                named tuples derived from *self.sequence_file.chromsizes*.
                If *None*, estimates itself, may affect performance.
            verbose (bool, default False): if True, prints tqdm progressbar
                of processing granges in split mode.

        Returns:
            None

        Raises:
            ValueError: An error in cases the extraction mode is not
                one of *'bulk'*, *'split'*.
        """
        if not isinstance(outdir, SequencePath):
            outdir = SequencePath(outdir)
        if mode == 'split':
            with self.sequence_file:
                filenames = list()
                for grange in tqdm(self,
                                   desc='Getting fasta',
                                   unit='grange',
                                   disable=not verbose):
                    grange.sequence_file_path = outdir / (str(outfileprefix) +
                                                          grange.sequence_header +
                                                          '.fasta')
                    with open(grange.sequence_file_path, 'w') as output:
                        self.sequence_file.get_fasta_by_coord(grange,
                                                              output,
                                                              chromsizes=chromsizes)
                    filenames.append(grange.sequence_file_path)
        elif mode == 'bulk':
            with self.sequence_file:
                outfilename = outdir / (str(outfileprefix) + '.fasta')
                with open(outfilename, 'w') as output:
                    for grange in self:
                        self.sequence_file.get_fasta_by_coord(grange,
                                                              output,
                                                              chromsizes=chromsizes)
                filenames = [outfilename]
        else:
            raise ValueError(f"get_fasta mode {mode} is not "
                             f"one of ['split', 'bulk'].")
        return filenames

    def relation_mapping(self, other, mapping, relation, verbose=False):
        """Maps relations of two genomic ranges lists.

        Provided with the map of genomic ranges names in
        the form of the dict with keys as names of genomic
        ranges in *self* and values as lists of names of
        genomic ranges in *other*, fills a *relation* for each
        genomic range in *self*, that's name is in the map keys,
        with a list of corresponding genomic ranges in other.

        Args:
            other (GenomicRangesList): the other genomic ranges
                list to fill relations of self genomic ranges list.
            mapping (dict of lists): mapping of genomic ranges names
                in *self* to the list of genomic ranges names in *other*.
            relation (str): the relation to fill.
            verbose (bool, default False): whether to report progress
                via tqdm progress bar or not.

        Returns:
            None
        """
        for grange in tqdm(self,
                           desc='Mapping relations',
                           unit='grange',
                           disable=not verbose):
            relating_granges = list()
            if grange.name in mapping.keys():
                values = mapping[grange.name]
                for code in values:
                    try:
                        other_granges = other.name_mapping[code]
                        relating_granges += [i for i in other_granges]
                    except KeyError:
                        continue
            if grange.relations.get(relation) is None:
                grange.relations[relation] = GenomicRangesList(relating_granges,
                                                               other.sequence_file_path)
            else:
                grange.relations[relation].update(relating_granges)

    @classmethod
    def parse_annotation(cls,
                         fileobj,
                         fileformat=None,
                         sequence_file_path=None,
                         column_names=None,
                         dtypes=None,
                         name_column=None,
                         name_pattern=None,
                         comment=None,
                         sep=None,
                         start_type=None,
                         end_type=None):
        """Parses annotation file into GenomicRangesList.

        The class method for loading GenomicRangesList instance
        from genomic annotation file in gtf, gff, bed3, bed6
        and bed12 formats. User can provide one's own file
        formats with given parsing patterns: column_names,
        dtypes, name_pattern, name_column, comment and sep. These
        parsing patterns are already provided for expected file
        formats. By default they are all None, so the class automatically
        chooses them by provided fileformat. If fileformat is None, expects
        one to provide all parsing patterns. Expects file with no header,
        but possibly with commented lines.

        Args:
            fileobj (fileobj): File-like object to parse
                annotation from.
            fileformat (str, default None): Type of file. Preferrably one
                of 'gtf', 'gff', 'bed3', 'bed6', 'bed12'. If
                None, means custom fileformat, so column_name,
                dtypes, name_pattern, comment and sep arguments must
                be provided.
            sequence_file_path (str, Path, SequencePath, default None):
                A path to the corresponding genome file.
            column_names (list of str, default None): List of column names. Has
                to contain 'chrom', 'start', 'end' values.
            dtypes (list of callables, default None): List of types of column values.
                Each callable corresponds to one column.
            name_column (str, default None): Name of column containing the name
                of each genomic range. May be not specified so
                the genomic range name will be derived from its coordinates.
            name_pattern (r-str, default None): Regexp pattern to extract genomic range
                name from column specified in name_column. Nay be not
                specified so the genomic range name will be derived from
                its coordinates.
            comment (str, default None): A character being used to comment lines.
            sep (str, default None): A character used to separate fields in the
                record (i.e. columns).
            start_type (int, default None): Chromosome coordinate start type:
                0-based or 1-based. One of: 0, 1.
            end_type (str, default None): Chromosome coordinate end type: inclusive or
                exclusive. Inclusive means character under the end coordinate
                includes in the sequence, exclusive means the opposite. One
                of: 'exclusive', 'inclusive'.

        Returns:
            GenomicRangesList: A list of genomic ranges corresponding
            to the given genome annotation.

        Raises:
            ValueError: An error if fileformat not one of supported (i.e
                'gtf', 'gff', 'bed3', 'bed6', 'bed12', None).
            ValueError: An error if fileformat is None and one of annotation
                patterns is None as well.
        """
        if fileformat not in cls.fileformats:
            raise ValueError(f"fileformat {fileformat} not one "
                             f"of {cls.fileformats}.")
        parser_dict = {'column_names': column_names,
                       'dtypes': dtypes,
                       'name_pattern': name_pattern,
                       'name_column': name_column,
                       'comment': comment,
                       'sep': sep,
                       'start_type': start_type,
                       'end_type': end_type}
        class_attr_map = {'column_names': 'column_names',
                          'dtypes': 'dtypes',
                          'name_pattern': 'name_patterns',
                          'name_column': 'name_columns',
                          'comment': 'comments',
                          'sep': 'seps',
                          'start_type': 'start_types',
                          'end_type': 'end_types'}
        for pattern_name, pattern in parser_dict.items():
            if pattern is None:
                if fileformat is None:
                    raise ValueError(f"Custom fileformat provided as None "
                                     f"and no custom {pattern_name} "
                                     f"provided as well.")
                parser_dict[pattern_name] = cls.__dict__[class_attr_map[pattern_name]][fileformat]
        if parser_dict['end_type'] not in ('inclusive', 'exclusive'):
            raise ValueError(f"end_type argument `{parser_dict['end_type']}` "
                             f"is not one of ['inclusive', 'exclusive'].")
        if parser_dict['start_type'] not in (0, 1):
            raise ValueError(f"start_type argument `{parser_dict['start_type']}' "
                             f"is not one of [0, 1].")
        record_holder = ([dtype(item)
                          for dtype, item in zip(parser_dict['dtypes'],
                                                 line.strip().split(parser_dict['sep']))]
                         for line in fileobj
                         if not line.startswith(parser_dict['comment']))
        annotated_holder = (dict(zip(parser_dict['column_names'],
                                     record))
                            for record in record_holder)
        grange_list = list()
        for record in annotated_holder:
            name = re.search(parser_dict['name_pattern'],
                             record.get('data', ""))
            name = name.group(1) if name is not None else None
            if parser_dict['start_type'] == 1:
                record['start'] -= 1
                record['end'] -= 1
            if parser_dict['end_type'] == 'inclusive':
                record['end'] += 1
            grange_list.append(GenomicRange(name=name,
                                            genome=sequence_file_path,
                                            **record))
        return cls(collection=grange_list,
                   sequence_file_path=sequence_file_path)


def extract_taxid_mapping(mapping, taxid):
    """Extracts mapping of orthologs for one taxid.

    Given the dictionary produced by `ortho2align
    get_orthodb_map` retrieves mapping of orthologs
    from query species to the provided taxid.

    Args:
        mapping (dict): mapping of orthologs
            loaded from json file produced
            by `ortho2align get_orthodb_map`.
        taxid (int, str): NCBI taxid of subject
            species.

    Returns:
        (dict) dictionary with query species
            gene ids as keys and list of subject
            species orthologs' ids as values.
    """
    taxid_mapping = dict()
    taxid = str(taxid)
    for key, value in mapping.items():
        try:
            taxid_mapping[key] = value[taxid]
        except KeyError:
            continue
    return taxid_mapping
