import re
from pathlib import Path
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import namedtuple, defaultdict
from sortedcontainers import SortedKeyList
from .alignment_utils import Alignment
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


class LocationNotSpecified(GenomicException):
    """Occurs when location of sequence file is not specified, but is accessed."""

    def __init__(self, instance):
        self.instance = instance
        super().__init__(f"Sequence file location for {self.instance} "
                         f"is not specified.")


def fasta_reformatter(infile, outfile, total, linewidth=60):
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

    Returns:
        None
    """
    read_counter = 0
    line_break = False
    while read_counter < total:
        c = infile.read(1)
        if c == "\n":
            continue
        else:
            read_counter += 1
            outfile.write(c)
            line_break = False
            if read_counter % linewidth == 0:
                outfile.write('\n')
                line_break = True
    if not line_break:
        outfile.write('\n')


class GenomicRange:
    """Represents one genomic range.

    Attributes:
        chrom
        start
        end
        strand
        _sequence_file_loc
        sequence_header
        name
        relations
    """

    def __init__(self, chrom, start, end, strand='.', name=None,
                 genome=None, sequence_file_loc=None, relations=dict(),
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
            sequence_file_loc: location of genomic range
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
        self._sequence_file_loc = sequence_file_loc
        self.sequence_header = f"{self.chrom}:{self.start}" \
                               f"-{self.end}{self.strand}"
        self.name = self.sequence_header if name is None else name
        self.relations = relations

    def __repr__(self):
        """GenomicRange instance representation."""
        return f"GenomicRange({self.chrom}, {self.start}, " \
               f"{self.end}, strand={self.strand}, name={self.name}, " \
               f"genome={self.genome}, sequence_file_loc=" \
               f"{self._sequence_file_loc}, relations={self.relations})"

    def __str__(self):
        """GenomicRange instance string representation."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}\t" \
               f"{self.name}\t{self.genome}"

    def to_dict(self):
        """Returns dict representation of the genomic range."""
        return {'chrom': self.chrom,
                'start': self.start,
                'end': self.end,
                'strand': self.strand,
                'genome': self.genome,
                'sequence_file_loc': self._sequence_file_loc,
                'name': self.name,
                'relations': {key: value.to_dict()
                              for key, value in self.relations}}

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
                         ('sequence_file_loc', None),
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

    def align(self, other, **kwargs):
        """Aligns two genomic ranges with blastn.

        If both genomic ranges have their sequences
        stored in separate files, aligns them with
        standalone blastn and specified parameters.
        Alignment file format is "7 std score".
        Standalone blast package should be installed
        manually.

        Args:
            other (GenomicRange): genomic range to
                align self with.
            kwargs (dict): any command line arguments
                to blastn in the following scheme:
                {'flag': 'value'}.

        Returns:
            (AlignedRangePair): aligned range pair,
                which associates self and other
                genomic ranges and their alignment.

        Raises:
            LocationNotSpecified in case one or both
            genomic ranges do not have separate
            sequence file.
        """
        if self._sequence_file_loc is None:
            raise LocationNotSpecified(self)
        if other._sequence_file_loc is None:
            raise LocationNotSpecified(other)
        command = ['blastn',
                   '-task', 'blastn',
                   '-query', self._sequence_file_loc,
                   '-subject', other._sequence_file_loc,
                   '-outfmt', '"7 std score"']
        add_args = sum((['-' + key, value]
                        for key, value in kwargs.items()),
                       [])
        with Popen(command + add_args, stdout=PIPE) as proc:
            alignment = Alignment.from_file(TextIOWrapper(proc.stdout))
        return AlignedRangePair(self, other, alignment)

    def align_with_relations(self, relation='synteny', **kwargs):
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
            alignment_list.append(self.align(grange, **kwargs))
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
        return f"AlignedRangePair({repr(self.query_range)}, " \
               f"{repr(self.subject_range)}, {repr(self.alignment)})"

    def __str__(self):
        """Returns string representation."""
        return f"AlignedRangePair\n{self.query_grange}\n" \
               f"{self.subject_grange}\n{self.alignment}"

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
        location (Path): path to the file.
        chromsizes (dict): dictionary of
            ChromosomseLocation instances.
        _file_obj: opened file under the
            location.
    """

    def __init__(self, location):
        """Initializes FastaSeqFile instance.

        Args:
            location (str, Path): path to the file.
        """
        self.location = Path(location)
        self._chromsizes = dict()
        self._file_obj = None

    def __repr__(self):
        """Returns representation of FastaSeqFile object."""
        return f"FastaSeqFile({self.location})"

    def __str__(self):
        """Returns string representation."""
        return f"FastaSeqFile({self.location})"

    def to_dict(self):
        """Returns dict representation."""
        return {'location': str(self.location)}

    @classmethod
    def from_dict(cls, dict_):
        """Restores FastaSeqFile instance from its dict representation."""
        keys = ['location']
        try:
            kwargs = {key: dict_[key] for key in keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of nececcary keys {keys}.")
        return cls(**kwargs)

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
            LocationNotSpecified in case the location of
                sequence file is not specified.
        """
        if not self._chromsizes:
            if self.location is None:
                raise LocationNotSpecified(self)
            with open(self.location, 'r') as infile:
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

    def __enter__(self):
        """Context manager for opening location file."""
        if self.location is None:
            raise LocationNotSpecified(self)
        self._file_obj = open(self.location, 'r')

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager for closing location file."""
        self._file_obj.close()

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
                                          self.location,
                                          self.chromsizes.keys())
        if grange.end > chrom_size:
            raise GenomicCoordinatesError(grange,
                                          self.chromsizes[grange.chrom])
        self._file_obj.seek(chrom_start + grange.start, 0)
        outfile.write(">" + grange.sequence_header + "\n")
        fasta_reformatter(self._file_obj,
                          outfile,
                          total=(grange.end - grange.start),
                          linewidth=linewidth)


class GenomicRangesList(SortedKeyList):
    """Represents a list of genomic ranges.

    Attributes:
        genome (str, Path): path to the corresponding genome file.
        source (FastaSeqFile): FastaSeqFile instance with the
            corresponding genome file.
        name_mapping (dict of lists): dictionary with keys
            as genomic ranges names and values as lists of
            corresponding GenomicRange instances.

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
        dtypes (dict of clists): dict of lists of callables representing
            dtypes of columns for selected file formats of genomic
            annotations.
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

    def __init__(self, collection=[], genome=None):
        """Initializes GenomicRangesList instance.

        Args:
            collection (iterable): collection of GenomicRange
                instances with the same genome.
            genome (str, Path): path to corresponding genome
                file.
        Returns:
            None.
        """
        super().__init__(iterable=collection,
                         key=lambda x: (x.chrom, x.start, x.end))
        self.genome = genome
        self.source = FastaSeqFile(self.genome)
        self._name_mapping = defaultdict(list)

    def __repr__(self):
        """Returns representaion."""
        granges = ", ".join(repr(grange) for grange in self)
        return f"GenomicRangesList([{granges}], genome={self.genome})"

    def __str__(self):
        """Returns string representation."""
        granges = "\n".join(str(grange) for grange in self)
        return f"GenomicRangesList of {self.genome}\n" \
               f"{granges}"

    def to_dict(self):
        """Returns dict representation."""
        return {'collection': [item.to_dict() for item in self],
                'genome': str(self.genome)}

    @classmethod
    def from_dict(cls, dict_):
        """Restores GenomicRangesList instance from dict representation."""
        keys = ['collection',
                'genome']
        try:
            kwargs = {key: dict_[key] for key in keys}
        except KeyError as e:
            raise ValueError(f"Key {e} is not in the list "
                             f"of nececcary keys {keys}.")
        kwargs['collection'] = [GenomicRange.from_dict(item)
                                for item in kwargs['collection']]
        return cls(**kwargs)

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
        new_range_list = GenomicRangesList([], self.genome)
        if len(self) == 0:
            return new_range_list
        new_range_list.add(self[0])
        for grange in self[1:]:
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
        for i in range(len(merged) - 1):
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
        return GenomicRangesList(inverted_granges, self.genome)

    def get_neighbours(self, other, distance=0):
        """Finds closest genomic ranges of self in other at the given distance.

        Args:
            other (GenomicRangesList): the other genomic
                ranges list to find neighbours of self in.
            distance (int): the distance at which to consider
                genomic ranges to be neighbours.

        Returns:
            None

        Raises:
            InconsistentGenomesError in case genomic ranges
                have inconsistent genomes.
        """
        self_index, other_index = 0, 0
        current_self, current_other = self_index, other_index
        while self_index < len(self) and other_index < len(other):
            try:
                if self[current_self].distance(other[current_other]) < - distance:
                    current_other += 1
                elif abs(self[current_self].distance(other[current_other])) <= distance:
                    self[current_self].relations['neighbours'].add(other[current_other])
                    current_other += 1
                else:
                    self_index += 1
                    current_self = self_index
                    current_other = other_index
            except InconsistentChromosomesError:
                self_index = current_self
                other_index = current_other
                if self[self_index].chrom < other[other_index].chrom:
                    self_index += 1
                    current_self = self_index
                else:
                    other_index += 1
                    current_other = other_index
            except IndexError:
                if current_other >= len(other):
                    current_other = other_index
                    current_self += 1
                    self_index = current_self
                else:
                    self_index += 1
                    current_self = self_index
                    current_other = other_index

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
        with self.source:
            if mode == 'split':
                filenames = list()
                for grange in self:
                    grange._sequence_file_loc = (str(outfileprefix) +
                                                 grange.sequence_header +
                                                 '.fasta')
                    with open(grange._sequence_file_loc, 'w') as output:
                        self.source.get_fasta_by_coord(grange, output)
                    filenames.append(grange._sequence_file_loc)
                return filenames
            elif mode == 'bulk':
                outfilename = str(outfileprefix) + '.fasta'
                with open(outfilename, 'w') as output:
                    for grange in self:
                        self.source.get_fasta_by_coord(grange, output)
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
        for key, value in mapping.items():
            try:
                self_grange = self.name_mapping[key]
            except KeyError:
                continue
            if self_grange.relations.get(relation) is None:
                self_grange.relations[relation] = GenomicRangesList([],
                                                                    other.genome)
            for code in value:
                try:
                    other_grange = other.name_mapping[code]
                    self_grange.relations[relation].append(other_grange)
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
        with NonExceptionalProcessPool(cores,
                                       verbose,
                                       suppress_exceptions) as p:
            alignments, exceptions = p.map(lambda x: x.align_with_relations(relation,
                                                                            **kwargs),
                                           self)

        return alignments, exceptions

    @classmethod
    def parse_annotation(cls,
                         fileobj,
                         fileformat=None,
                         genome=None,
                         column_names=None,
                         dtypes=None,
                         name_column=None,
                         name_pattern=None,
                         comment=None,
                         sep=None):
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
            genome (str, Path): path to the corresponding genome file
                (default: None).
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
                       'sep': sep}
        for pattern_name, pattern in parser_dict.items():
            if pattern is None:
                if fileformat is None:
                    raise ValueError(f"Custom fileformat provided as None "
                                     f"and no custom {pattern_name} "
                                     f"provided as well.")
                parser_dict[pattern_name] = cls.__dict__[pattern_name + 's'][fileformat]
        record_holder = ([dtype(item)
                          for dtype, item in zip(parser_dict['dtypes'],
                                                 line.strip().split(parser_dict['sep']))]
                         for line in fileobj
                         if not line.startswith(parser_dict['comment']))
        annotated_holder = (dict(zip(parser_dict['column_names'],
                                     record))
                            for record in record_holder)
        grangeslist = GenomicRangesList([], genome=genome)
        for record in annotated_holder:
            name = re.search(parser_dict['name_pattern'],
                             record.get('data', "")).group(1)
            name = name if name else None
            grangeslist.add(GenomicRange(name=name, genome=genome, **record))
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
