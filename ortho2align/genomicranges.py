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
from functools import partial
from itertools import chain
from pathlib import Path
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import namedtuple, defaultdict
from sortedcontainers import SortedKeyList
from tqdm import tqdm
from .alignment import Alignment
from .parallel import NonExceptionalProcessPool


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


def fasta_reformatter(infile, outfile, total, linewidth=60,
                      reverse=False, complement=False):
    """Reformats fasta file.

    Allows one to reformat existing sequence in
    fasta format to the specified linewidth.

    Args:
        infile (fileobj): input fasta file with cursor set
            at the start of sequence fragment.
        outfile (fileobj): output fasta file.
        total (int): number of symbols to read.
        linewidth (int): linewidth of new sequence
            (default: 60).
        reverse (bool): if True, takes sequence backwards
            (default: False).
        complement (bool): if True, returns complement
            nucleotide sequence (default: False).

    Returns:
        None
    """
    read_counter = 0
    line_break = False
    compl_map = {'A': 'T',
                 'a': 't',
                 'T': 'A',
                 't': 'a',
                 'G': 'C',
                 'g': 'c',
                 'C': 'G',
                 'c': 'g'}
    while read_counter < total:
        c = infile.read(1)
        if reverse:
            infile.seek(infile.tell() - 2)
        if c == "\n":
            continue
        else:
            if complement:
                c = compl_map[c]
            read_counter += 1
            outfile.write(c)
            line_break = False
            if read_counter % linewidth == 0:
                outfile.write('\n')
                line_break = True
    if not line_break:
        outfile.write('\n')


class SequencePath:
    """Wrapper aroudn pathlib.Path with checking of None.

    Attributes:
        path (None or Path): filesystem path.
    """

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

    def to_dict(self):
        """Returns dict representation."""
        if self.path is None:
            return None
        else:
            return str(self.path)

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


class GenomicRange:
    """Represents one genomic range.

    Genomic coordinates are 0-based exclusive:
    * chromosome start coordinate is 0;
    * character under genomic range coordinate
        is not included into the sequence.

    Attributes:
        chrom (str): name of chromosome.
        start (int): start of genomic range.
        end (int): end of genomic range.
        strand (str): genomic range strand. One of
            "+", "-", "." (inessential strandness)
            (default: None).
        name (str): name of genomic range (defalut: None).
        genome (str): corresponding genome filename.
        sequence_file_path: location of genomic range
            sequence file.
        relations (dict): a dictionary with related GenomicRangesList
            instances.
    """

    def __init__(self, chrom, start, end, strand='.', name=None,
                 genome=None, sequence_file_path=None, relations=None,
                 **kwargs):
        """Initializes GenomicRange instance.

        Args:
            chrom (str): name of chromosome.
            start (int): start of genomic range.
            end (int): end of genomic range.
            strand (str): genomic range strand. One of
                "+", "-", "." (inessential strandness)
                (default: None).
            name (str): name of genomic range (defalut: None).
            genome (str): corresponding genome filename.
            sequence_file_path: location of genomic range
                sequence file.
            relations (dict): a dictionary with related GenomicRangesList
                instances.
            kwargs (dict): other arguments.

        Returns:
            None.
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.genome = genome
        self._sequence_file_path = SequencePath(sequence_file_path)
        self._name = None if name is None else name
        self.relations = dict() if relations is None else relations

    def __repr__(self):
        """GenomicRange instance representation."""
        return f"GenomicRange({self.chrom}, {self.start}, " \
               f"{self.end}, strand={self.strand}, name={self.name}, " \
               f"genome={self.genome}, sequence_file_path=" \
               f"{self.sequence_file_path}, relations={self.relations})"

    def __str__(self):
        """GenomicRange instance string representation."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}\t" \
               f"{self.name}\t{self.genome}"

    def __eq__(self, other):
        """Equality magic method."""
        return (self.chrom == other.chrom and self.start == other.start and
                self.end == other.end and self.strand == other.strand and
                self.genome == other.genome and
                self.sequence_file_path == other.sequence_file_path and
                self.name == other.name and self.relations == other.relations)

    def to_dict(self):
        """Returns dict representation of the genomic range."""
        return {'chrom': self.chrom,
                'start': self.start,
                'end': self.end,
                'strand': self.strand,
                'genome': self.genome,
                'sequence_file_path': self.sequence_file_path.to_dict(),
                'name': self.name,
                'relations': {key: value.to_dict()
                              for key, value in self.relations.items()}}

    @classmethod
    def from_dict(cls, dict_):
        """Restores genomic range from dict representation.

        Args:
            dict_: dict representation of GenomicRange
                instance generated with `GenomicRange.to_dict()`.

        Returns:
            (GenomicRange) restored GenomicRange instance.
        """
        arbitrary_keys = ['chrom',
                          'start',
                          'end']
        optional_keys = [('strand', '.'),
                         ('name', None),
                         ('genome', None),
                         ('sequence_file_path', None),
                         ('relations', dict())]
        try:
            kwargs = {key: dict_[key] for key in arbitrary_keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of arbitrary keys {arbitrary_keys}.")
        kwargs.update({key[0]: dict_.get(key[0], key[1])
                       for key in optional_keys})
        if kwargs['relations']:
            kwargs['relations'] = {key: GenomicRangesList.from_dict(value)
                                   for key, value in kwargs['relations'].items()}
        return cls(**kwargs)

    @property
    def sequence_file_path(self):
        """`sequence_file_path` getter."""
        return self._sequence_file_path

    @sequence_file_path.setter
    def sequence_file_path(self, value):
        """`sequence_file_path` setter."""
        self._sequence_file_path = SequencePath(value)

    @property
    def sequence_header(self):
        return f"{self.chrom}:{self.start}-{self.end}{self.strand}"

    @property
    def name(self):
        if self._name is None:
            return self.sequence_header
        return self._name

    def distance(self, other):
        """Calculates distance between two genomic ranges.

        In case of same genomes and chromosomes
        the distance between two genomic ranges
        can be calculated regardless of their
        strandness. If two genomic ranges overlap,
        the distance is zero. In case of inconsistent
        chromosomes or genomes the distance is undefined.

        Args:
            other (GenomicRange): the genomic range
                to which calculate the distance from
                self.

        Returns:
            (int) distance, non-negative number.

        Raises:
            InconsistentGenomesError in case of inconsistent
                genomes.
            InconsistentChromosomesError in case of inconsistent
                chromosomes.
        """
        if self.genome != other.genome:
            raise InconsistentGenomesError(self.genome, other.genome)
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        if self.end < other.start:
            return other.start - self.end
        if other.end < self.start:
            return other.end - self.start
        return 0

    def merge(self, other):
        """Merges two genomic ranges together.

        In case of the same genomes and chromosomes
        two genomic ranges can be merged strand-nonspecifically
        into one genomic range which starts where the closest to
        the chromosome start genomic range starts and ends where
        the closest to the chromosome end genomic range ends.

        Args:
            other (GenomicRange): genomic range which is merged
                with self.

        Returns:
            (GenomicRange) new genomic range with "." strand and
                the same genome.

        Raises:
            InconsistentGenomesError in case of inconsistent
                genomes.
            InconsistentChromosomesError in case of inconsistent
                chromosomes.
        """
        if self.genome != other.genome:
            raise InconsistentGenomesError(self.genome, other.genome)
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        return GenomicRange(chrom=self.chrom,
                            start=start,
                            end=end,
                            strand=".",
                            genome=self.genome)

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
            (GenomicRange) flanked genomic range.
        """
        if - flank_distance * 2 > self.end - self.start:
            raise ValueError(f"Flank distance value {flank_distance} "
                             "is negative and greater than the "
                             "half-length of the genomic range.")
        start = self.start - flank_distance
        end = self.end + flank_distance
        return GenomicRange(chrom=self.chrom,
                            start=start,
                            end=end,
                            strand=".",
                            genome=self.genome)

    def find_neighbours(self, other_list, distance=0):
        """Finds neighbours of self in other list at given distance.

        Args:
            other_list (GenomicRangesList): other list to find
                neighbours in.
            distance (int, float): distance at which two
                genomic ranges are considered neighbours
                (default: 0).

        Returns:
            (GenomicRangesList) genomic ranges list with all
            neighbours of self in other_list at given distance.
        """
        neighbours = GenomicRangesList([],
                                       sequence_file_path=other_list.sequence_file_path)
        index = other_list.bisect_left(self)
        left_index = index - 1
        right_index = index
        try:
            while abs(self.distance(other_list[left_index])) <= distance:
                neighbours.add(other_list[left_index])
                left_index -= 1
        except InconsistentChromosomesError:
            pass
        except IndexError:
            pass
        try:
            while abs(self.distance(other_list[right_index])) <= distance:
                neighbours.add(other_list[right_index])
                right_index += 1
        except InconsistentChromosomesError:
            pass
        except IndexError:
            pass
        return neighbours

    def align_blast(self, other, program='blastn', task='blastn', **kwargs):
        """Aligns two genomic ranges with BLAST.

        If both genomic ranges have their sequences
        stored in separate files, aligns them with
        standalone blastn and specified parameters.
        Alignment file format is "7 std score".
        Standalone blast package should be installed
        manually.

        Args:
            other (GenomicRange): genomic range to
                align self with.
            program (str): a program from Standalone
                BLAST to use (default: 'blastn')
            task (str): which task from a BLAST program
                use (default: 'blastn')
            kwargs (dict): any command line arguments
                to blastn in the following scheme:
                {'flag': 'value'}.

        Returns:
            (AlignedRangePair): aligned range pair,
                which associates self and other
                genomic ranges and their alignment.

        Raises:
            FilePathNotSpecifiedError in case one or both
            genomic ranges do not have separate
            sequence file.
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
        add_args = sum((['-' + key, value]
                        for key, value in kwargs.items()),
                       [])
        with Popen(command + add_args, stdout=PIPE, stderr=PIPE) as proc:
            alignment = Alignment.from_file_blast(TextIOWrapper(proc.stdout),
                                                  self.chrom,
                                                  other.chrom,
                                                  self.name,
                                                  other.name,
                                                  self.strand,
                                                  other.strand)
        return AlignedRangePair(self, other, alignment)

    def align_with_relation(self, relation, **kwargs):
        """Aligns genomic range with all its relations.

        Aligns given genomic range with all its
        relations using standalone blastn
        and returns list of AlignedRangePair
        instances.

        Args:
            relation (str): relation type.
            kwargs (dict): any command line arguments
                to blastn in the following scheme:
                {'flag': 'value'}.

        Returns:
            (list) list of AlignedRangePair instances.

        Raises:
            ValueError in case provided relation is not
                available for the genomic range.
        """
        if relation not in self.relations.keys():
            raise ValueError(f"Relation {relation} not in "
                             f"available list of relations: "
                             f"{self.relations.keys()}.")
        alignment_list = list()
        for grange in self.relations[relation]:
            alignment_list.append(self.align_blast(grange, **kwargs))
        return alignment_list

    def merge_relations(self, relation, distance=0):
        """Merges genomic ranges from relation on specified distance.

        All genomic ranges closer to each other then the
        specified distance are merged together and returned.

        Args:
            relation (str): relation type.
            distance (int): distance to merge
                genomic ranges in relation
                (default: 0).

        Returns:
            (GenomicRangesList) merged genomic ranges from
                relation.

        Raises:
            ValueError in case there is no given relation
                availabe for provided genomic range.
        """
        try:
            return self.relations[relation].merge(distance)
        except KeyError:
            raise ValueError(f"There is no {relation} "
                             f"relations for the {self}.")


def align_with_relation_wrapper(grange, relation, **kwargs):
    """Simple wrapper around GenomicRange method for multiprocessing."""
    return grange.align_with_relation(relation, **kwargs)


class AlignedRangePair:
    """An association of genomic ranges and their alignment.

    Stores query genomic range, subject genomic range and
    their alignment.

    Attributes:
        query_grange (GenomicRange): query genomic range.
        subject_grange (GenomicRange): subject genomic range.
        alignment (Alignment): alignment of query and subject
        genomic ranges.
    """

    def __init__(self, query_grange, subject_grange, alignment):
        """Initializes AlignedRangePair instance.

        Args:
            query_grange (GenomicRange): query genomic range.
            subject_grange (GenomicRange): subject genomic range.
            alignment (Alignment): alignment of query and subject
            genomic ranges.

        Returns:
            None.
        """
        self.query_grange = query_grange
        self.subject_grange = subject_grange
        self.alignment = alignment

    def __repr__(self):
        """Returns representation of AlignedRangePair instance."""
        return f"AlignedRangePair({repr(self.query_grange)}, " \
               f"{repr(self.subject_grange)}, {repr(self.alignment)})"

    def __str__(self):
        """Returns string representation."""
        return f"AlignedRangePair\n{self.query_grange}\n" \
               f"{self.subject_grange}\n{self.alignment}"

    def __eq__(self, other):
        """Equality magic method."""
        return (self.query_grange == other.query_grange and
                self.subject_grange == other.subject_grange and
                self.alignment == other.alignment)

    def to_dict(self):
        """Returns dict representation of AlignedRangePair instance."""
        return dict(query_grange=self.query_grange,
                    subject_grange=self.subject_grange,
                    alignment=self.alignment)

    @classmethod
    def from_dict(cls, dict_):
        """Restores AlignedRangePair instance from dict representation."""
        keys = ['query_grange',
                'subject_grange',
                'alignment']
        try:
            kwargs = {key: dict_[key] for key in keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of nececcary keys {keys}.")
        return cls(**kwargs)


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
        return f"FastaSeqFile({self.sequence_file_path})"

    def __str__(self):
        """Returns string representation."""
        return f"FastaSeqFile({self.sequence_file_path})"

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
            kwargs = {key: dict_[key] for key in keys}
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

    def locate_coordinate(self, chrom_start, coord):
        """Puts file pointer to the beginning of the genomic range.

        Args:
            chrom_start (int): byte number of chromosome start.
            coord (int): coordinate in the corresponding chromosome.

        Returns:
            None
        """
        self._file_obj.seek(chrom_start + coord + coord // self.line_length)

    def get_fasta_by_coord(self, grange, outfile, linewidth=60):
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
            linewidth (int): line length of fasta sequence to write
                (default: 60).

        Returns:
            None.

        Raises:
            ChromosomeNotFoundError if the chromosome of
                the given genomic range is not presented
                in FastaSeqFile instance.
        """
        try:
            chrom_size, chrom_start = self.chromsizes[grange.chrom]
        except KeyError:
            raise ChromosomeNotFoundError(grange.chrom,
                                          self.sequence_file_path,
                                          self.chromsizes.keys())
        if grange.end > chrom_size:
            raise GenomicCoordinatesError(grange,
                                          self.chromsizes[grange.chrom])
        if grange.strand == '-':
            self.locate_coordinate(chrom_start, grange.end - 1)
            reverse = True
            complement = True
        else:
            self.locate_coordinate(chrom_start, grange.start)
            reverse = False
            complement = False
        outfile.write(">" + grange.sequence_header + "\n")
        fasta_reformatter(self._file_obj,
                          outfile,
                          total=(grange.end - grange.start),
                          linewidth=linewidth,
                          reverse=reverse,
                          complement=complement)


class GenomicRangesList(SortedKeyList):
    """Represents a list of genomic ranges.

    Genomic coordinates are 0-based exclusive:
    * chromosome start coordinate is 0;
    * character under genomic range coordinate
        is not included into the sequence.

    Attributes:
        sequence_file_path (str, Path): path to the corresponding
            genome file.
        sequence_file (FastaSeqFile): FastaSeqFile instance with the
            corresponding genome file.
        name_mapping (dict of lists): dictionary with keys
            as genomic ranges names and values as lists of
            corresponding GenomicRange instances.
        _verbose_amount (int): how many genomic ranges to
            print at maximum (positive number) (default: 10).

    Class attributes:
        column_names (dict of lists): column names of selected
            file formats of genomic annotations.
        name_columns (dict of strs): names of columns containing
            information about genomic range names for selected
            file formats of genomic annotations.
        comments (dict of strs): comment characters for selected file
            formats of genomic annotations.
        name_patterns (dict of r-strs): regex patterns to extract
            genomic range from the column containing the name for
            selected file formats of genomic annotations.
        seps (dict of strs): field separator characters for selected
            file formats of genomic annotations.
        fileformats (list): list of selected file formats of genomic
            annotations (gtf, gff, bed3, bed6, bed12, None -- user-specified).
        dtypes (dict of lists): dict of lists of callables representing
            dtypes of columns for selected file formats of genomic
            annotations.
        start_types (dict of ints): maps file format to genomic range
            coordinate start type (0 or 1).
        end_types (dict of strs): maps file format to genomic range
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
                             'name',
                             'score',
                             'strand'],
                    'bed12': ['chrom',
                              'start',
                              'end',
                              'name',
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
                    'bed6': 'name',
                    'bed12': 'name'}
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
                       int,
                       str],
              'bed12': [str,
                        int,
                        int,
                        str,
                        int,
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

    @classmethod
    def key_func(cls, grange):
        """Returns key for sorting of genomic ranges."""
        return (grange.chrom, grange.start, grange.end)

    def __new__(cls, collection=[], sequence_file_path=None):
        """Makes new instance of GenomicRangesList."""
        instance = super().__new__(cls,
                                   iterable=collection,
                                   key=cls.key_func)
        return instance

    def __init__(self, collection=[], sequence_file_path=None):
        """Initializes GenomicRangesList instance.

        Args:
            collection (iterable): collection of GenomicRange
                instances with the same genome.
            sequence_file_path (str, Path): path to corresponding genome
                file.
        Returns:
            None.
        """
        super().__init__(iterable=collection,
                         key=GenomicRangesList.key_func)
        self._sequence_file_path = SequencePath(sequence_file_path)
        self.sequence_file = FastaSeqFile(self.sequence_file_path)
        self._name_mapping = defaultdict(list)
        self._verbose_amount = 10

    def __repr__(self):
        """Returns representaion."""
        granges = ", ".join([repr(grange) for grange in self][:self._verbose_amount])
        appendix = ", ..." if self._verbose_amount < len(self) else ""
        return f"GenomicRangesList([{granges}{appendix}], " \
               f"sequence_file_path={self.sequence_file_path})"

    def __str__(self):
        """Returns string representation."""
        granges = "\n".join([str(grange) for grange in self][:self._verbose_amount])
        appendix = "\n..." if self._verbose_amount < len(self) else ""
        return f"GenomicRangesList of {len(self)} genomic ranges of"  \
               f"{self.sequence_file_path}\n{granges}{appendix}\n" \
               f"{min(self._verbose_amount, len(self))} of {len(self)} genomic ranges."

    def __eq__(self, other):
        """Equality magic method."""
        if len(self) != len(other):
            return False
        if self.sequence_file_path != other.sequence_file_path:
            return False
        return all([self_grange == other_grange
                    for self_grange, other_grange in zip(self, other)])

    def to_dict(self):
        """Returns dict representation."""
        return {'collection': [item.to_dict() for item in self],
                'sequence_file_path': str(self.sequence_file_path)}

    @classmethod
    def from_dict(cls, dict_):
        """Restores GenomicRangesList instance from dict representation."""
        keys = ['collection',
                'sequence_file_path']
        try:
            kwargs = {key: dict_[key] for key in keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of nececcary keys {keys}.")
        kwargs['collection'] = [GenomicRange.from_dict(item)
                                for item in kwargs['collection']]
        return cls(**kwargs)

    @property
    def sequence_file_path(self):
        """sequence_file_path getter."""
        return self._sequence_file_path

    @sequence_file_path.setter
    def sequence_file_path(self, value):
        """sequence_file_path setter."""
        self._sequence_file_path = SequencePath(value)

    @property
    def name_mapping(self):
        """GenomicRange name mapping.

        Dictionary that maps names of GenomicRange
            instances to the instances themselves.
        """
        if not self._name_mapping:
            for grange in self:
                self._name_mapping[grange.name].append(grange)
        return self._name_mapping

    def merge(self, distance=0):
        """Merges genomic ranges on the given distance.

        If the distance between two genomic ranges
        is smaller then the provided value, they are
        merged via GenomicRange.merge and added to
        the new GenomicRangesList instance.

        Args:
            distance: distance to merge proximal
                genomic ranges (default: 0).

        Returns:
            (GenomicRangeList) new genomic ranges list
                with all proximal genomic ranges merged.
                If the initial genomic ranges list is empty,
                returns empty genomic ranges list.
        """
        new_range_list = GenomicRangesList([], self.sequence_file_path)
        if len(self) == 0:
            return new_range_list
        new_range_list.add(self[0])
        for grange in tqdm(self[1:]):
            try:
                if abs(new_range_list[-1].distance(grange)) <= distance:
                    old_range = new_range_list.pop()
                    new_range_list.add(old_range.merge(grange))
                else:
                    new_range_list.add(grange)
            except InconsistentChromosomesError:
                new_range_list.add(grange)
        return new_range_list

    def inter_ranges(self, distance=0):
        """Returns inner genomic regions covered with no genomic ranges.

        First merges all genomic ranges provided
        in self at the given distance with
        `GenomicRangesList.merge`, then inverses
        them so as new regions do not contain any
        of the previous genomic ranges. Regions
        at the beginning and the end of chromosomes
        are not included.

        Args:
            distance: distance to merge initial
                genomic ranges (default: 0).

        Returns:
            (GenomicRangesList): new genomic ranges
                list with inverted genomic ranges.
        """
        merged = self.merge(distance)
        inverted_granges = list()
        for i in tqdm(range(len(merged) - 1)):
            if merged[i].chrom == merged[i + 1].chrom:
                new_chrom = merged[i].chrom
                new_start = merged[i].end
                new_end = merged[i + 1].start
                new_strand = "."
                new_genome = merged[i].genome
                inverted_granges.append(GenomicRange(new_chrom,
                                                     new_start,
                                                     new_end,
                                                     new_strand,
                                                     genome=new_genome))
        return GenomicRangesList(inverted_granges, self.sequence_file_path)

    def get_neighbours(self, other, distance=0, relation='neighbours'):
        """Finds all neighbours of self in other at the given distance.

        For each genomic range in self finds all ranges in other
        that are closer to the range than provided distance,
        forms new GenomicRangesList instance with them and assigns
        it to the `.relation` dict under the provided relation key.

        Args:
            other (GenomicRangesList): the other genomic ranges
                list in which to find neighbours of self.
            distance (int, float): distance at which two genomic
                ranges are considered neighbours (default: 0).
            relation (str): name of relation to assign all
                found neighbours of each genomic range in self
                (default: 'neighbours').

        Returns:
            None.
        """
        for grange in tqdm(self):
            grange.relations[relation] = grange.find_neighbours(other,
                                                                distance=distance)

    def flank(self, distance=0, check_boundaries=True):
        """Flanks all genomic ranges in the list with given distance.

        Allows correction of new start and end coordinates in case
        they exceeds corresponding chromosome size.

        Args:
            distance (int): value to flank genomic ranges
                (default: 0).
            check_boundaries (bool): whether to check if
                new start coordinates are not less than zero
                and new end coordinates are not greater than
                corresponding chromosome sizes. If True, corrects
                coordinates (default: False).

        Returns:
            (GenomicRangesList) new genomic ranges list with flanked
                and corrected genomic ranges.
        """
        new_grangelist = GenomicRangesList([],
                                           sequence_file_path=self.sequence_file_path)
        for grange in tqdm(self):
            new_grange = grange.flank(distance)
            if check_boundaries:
                new_grange.start = max(0, new_grange.start)
                new_grange.end = min(self.sequence_file.chromsizes[new_grange.chrom].size,
                                     new_grange.end)
            new_grangelist.add(new_grange)
        return new_grangelist

    def get_fasta(self, outfileprefix, mode='split'):
        """Extracts sequences of genomic ranges.

        Can extract sequences of all genomic ranges
        into single file named `outfileprefix` + '.fasta'
        or each sequence into separate files named
        `outfileprefix` + genomic_range.name + '.fasta'.

        Args:
            outfileprefix (str, Path): prefix of the
                output file(s).
            mode (str): the mode to extract sequences,
                one of 'split', 'bulk' (default: 'split').

        Returns:
            None

        Raises:
            ValueError in cases the extraction mode is not
                one of 'bulk', 'split'.
        """
        with self.sequence_file:
            if mode == 'split':
                filenames = list()
                for grange in tqdm(self):
                    grange.sequence_file_path = (str(outfileprefix) +
                                                 grange.sequence_header +
                                                 '.fasta')
                    with open(grange.sequence_file_path, 'w') as output:
                        self.sequence_file.get_fasta_by_coord(grange, output)
                    filenames.append(grange.sequence_file_path)
                return filenames
            elif mode == 'bulk':
                outfilename = str(outfileprefix) + '.fasta'
                with open(outfilename, 'w') as output:
                    for grange in self:
                        self.sequence_file.get_fasta_by_coord(grange, output)
                return [outfilename]
            else:
                raise ValueError(f"get_fasta mode {mode} not "
                                 f"one of ['split', 'bulk'].")

    def relation_mapping(self, other, mapping, relation):
        """Maps relations of two genomic ranges lists.

        Provided with the map of genomic ranges names in
        the form of the dict with keys as names of genomic
        ranges in self and values as lists of names of
        genomic ranges in other, fills a relation for each
        genomic range in self, that's name is in the map keys,
        with corresponding genomic ranges in other.

        Args:
            other (GenomicRangesList): the other genomic ranges
                list to fill relations of self genomic ranges list.
            mapping (dict of lists): mapping of genomic ranges names
                in self to list of genomic ranges names in other.
            relation (str): the relation to fill.

        Returns:
            None
        """
        for grange in tqdm(self):
            if grange.relations.get(relation) is None:
                grange.relations[relation] = GenomicRangesList([],
                                                               other.sequence_file_path)
            if grange.name in mapping.keys():
                values = mapping[grange.name]
                for code in values:
                    try:
                        other_granges = other.name_mapping[code]
                        grange.relations[relation].update(other_granges)
                    except KeyError:
                        continue

    def align_with_relation(self, relation, cores=1, verbose=False,
                            suppress_exceptions=False,
                            **kwargs):
        """Aligns each genomic range with its relation if it exists.

        For each genomic range in self in case its' relation
        exists, aligns that genomic range with each range in
        that relation using standalone blastn. Allows
        multiprocessing of the alignment process.

        Args:
            relation (str): relation to align with.
            cores (int): how many cores to use (default: 1).
            verbose (bool): if True, will print progress
                with tqdm progress bar (default: False).
            suppress_exceptions (bool): if True, will suppress
                reporting exceptions to the return list (dfault: False).
            kwargs (dict-like): any kwargs to pass to
                GenomicRange.align representing CLI arguments
                for blastn.
        Returns:
            (list, list) two lists: the first contains AlignedRangePair
                instances representing alignments, the cecond contains
                exceptions occured during alignment if `suppress_exceptions`
                is False (otherwise it is empty).

        Raises:
            EmptyGenomicRangesListError in case self genomic ranges
                list is empty.
        """
        if len(self) == 0:
            raise EmptyGenomicRangesListError(self)

        process_func = partial(align_with_relation_wrapper,
                               relation=relation,
                               **kwargs)

        with NonExceptionalProcessPool(cores,
                                       verbose,
                                       suppress_exceptions) as p:
            alignments, exceptions = p.map(process_func, self)

        if len(exceptions) > 0:
            print(exceptions)

        return list(chain.from_iterable(alignments))

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
            fileobj (fileobj): file-like object to parse
                annotation from.
            fileformat (str): type of file. Preferrably one
                of 'gtf', 'gff', 'bed3', 'bed6', 'bed12'. If
                None, means custom fileformat, so column_name,
                dtypes, name_pattern, comment and sep arguments must
                be provided (default: None).
            sequence_file_path (str, Path): path to the corresponding
                genome file (default: None).
            column_names (list of str): list of column names. Has
                to contain 'chrom', 'start', 'end' values (default: None).
            dtypes (list of callables): list of types of column values.
                Each callable corresponds to one column (default: None).
            name_column (str): name of column containing the name
                of each genomic range. May be not specified so
                the genomic range name will be derived from its coordinates
                (defalut: None).
            name_pattern (r-str): regex pattern to extract genomic range
                name from column specified in name_column. Nay be not
                specified so the genomic range name will be derived from
                its coordinates (default: None).
            comment (str): a character being used to comment lines
                (default: None).
            sep (str): a character used to separate fields in the
                record (i.e. columns) (default: None).
            start_type (int): chromosome coordinate start type: 0-based
                or 1-based. One of: 0, 1 (default: None).
            end_type (str): chromosome coordinate end type: inclusive or
                exclusive. Inclusive means character under the end coordinate
                includes in the sequence, exclusive means the opposite. One
                of: 'exclusive', 'inclusive' (default: None).

        Returns:
            (GenomicRangesList) list of genomic ranges corresponding
                to the given genome annotation.

        Raises:
            ValueError if fileformat not one of supported (i.e
                'gtf', 'gff', 'bed3', 'bed6', 'bed12', None).
            ValueError if fileformat is None and one of annotation
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
        grangeslist = cls(collection=[],
                          sequence_file_path=sequence_file_path)
        for record in annotated_holder:
            name = re.search(parser_dict['name_pattern'],
                             record.get('data', ""))
            name = name.group(1) if name is not None else None
            if parser_dict['start_type'] == 1:
                record['start'] -= 1
                record['end'] -= 1
            if parser_dict['end_type'] == 'inclusive':
                record['end'] += 1
            grangeslist.add(GenomicRange(name=name,
                                         genome=sequence_file_path,
                                         **record))
        return grangeslist


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
