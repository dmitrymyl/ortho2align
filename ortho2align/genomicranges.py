import json
import re
from pathlib import Path
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import namedtuple, defaultdict
from multiprocessing.pool import Pool
from sortedcontainers import SortedKeyList
from .alignment_utils import Alignment


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
        super().__init__(f"Empty GenomicRangesList instance: {self.range_list}")


class GenomicCoordinatesError(GenomicException):
    """Occurs when coordinates of GenomicRange instance are out of the chromosome."""

    def __init__(self, grange, chromosome):
        self.grange = grange
        self.chromosome = chromosome
        super().__init__(f"Coordinates of GenomicRange {self.grange} "
                         f"are out of chromosome with size {self.chromosome.size}")


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
                 genome=None, sequence_file_loc=None,
                 synteny=None, neighbours=None, **kwargs):
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
            synteny (GenomicRangesList): list of syntenic
                genomic ranges (default: None).
            neighbours (GenomicRangesList): list of neighbouring
                genomic ranges (default: None).
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
        # TODO: refactor relations. They should not be done this way
        # as we loses genome information.
        self.relations = {'synteny': (GenomicRangesList([], self.genome)
                                        if synteny is None else synteny),
                            'neighbours': (GenomicRangesList([], self.genome)
                                           if neighbours is None else neighbours)}

    def __repr__(self):
        """GenomicRange instance representation."""
        return f"GenomicRange({self.chrom}, {self.start}, "
               f"{self.end}, {self.strand})"

    def __str__(self):
        """GenomicRange instance string representation."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}"

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
            raise ValueError(f"relation {relation} not in "
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

    def to_dict(self):
        """Returns dict representation of the genomic range."""
        raise NotImplementedError

    @classmethod
    def from_dict(self, dict_):
        """Restores genomic range from dict representation."""
        raise NotImplementedError


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

    @property
    def chromsizes(self):
        """Poperty representing genomic chromosomes.

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
        genome
        source
        name_mapping

    Class attributes:
        annotation_patterns
        comments
        name_patterns
        seps
        filetypes
        dtypes
    """
    annotation_patterns = {'gff': ['chrom',
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
    filetypes = ['gtf', 'gff', 'bed3', 'bed6', 'bed12']
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
        super().__init__(iterable=collection,
                         key=lambda x: (x.chrom, x.start, x.end))
        self.genome = genome
        self.source = FastaSeqFile(self.genome)
        self._name_mapping = defaultdict(list)

    @property
    def name_mapping(self):
        if not self._name_mapping:
            for grange in self:
                self._name_mapping[grange.name].append(grange)
        return self._name_mapping

    def merge(self, distance=0):
        if len(self) == 0:
            raise EmptyGenomicRangesListError(self)
        new_range_list = GenomicRangesList([], self.genome)
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
                raise ValueError(f"get_fasta mode {mode} not one of ['split', 'bulk'].")

    def relation_mapping(self, other, mapping, relation):
        for key, value in mapping.items():
            for code in value:
                try:
                    self_grange = self.name_mapping[key]
                    other_grange = other.name_mapping[code]
                    self_grange.relations[relation].append(other_grange)
                except KeyError:
                    continue

    def align_with_relations(self, relation, cores=1):
        if len(self) == 0:
            raise EmptyGenomicRangesListError(self)
        if relation not in self[0].relations.keys():
            raise ValueError(f"relation type {relation} "
                             f"not in the list of available "
                             f"relations: {self[0].relations.keys()}.")
        with Pool(cores) as p:
            alignments = p.map(lambda x: x.align_with_relations(relation),
                               self)
        return alignments

    @classmethod
    def parse_annotation(cls,
                         fileobj,
                         filetype,
                         genome=None,
                         annotation_pattern=None,
                         dtypes=None,
                         name_pattern=None,
                         comment=None,
                         sep=None):
        if filetype not in cls.filetypes:
            raise ValueError(f"filetype {filetype} not one of {cls.filetypes}.")
        if annotation_pattern is None:
            annotation_pattern = cls.annotation_patterns[filetype]
        if dtypes is None:
            dtypes = cls.dtypes[filetype]
        if name_pattern is None:
            name_pattern = cls.name_patterns[filetype]
        if comment is None:
            comment = cls.comments[filetype]
        if sep is None:
            sep = cls.seps[filetype]
        record_holder = ([func(item)
                          for func, item in zip(dtypes,
                                                line.strip().split(sep))]
                         for line in fileobj
                         if not line.startswith(comment))
        annotated_holder = (dict(zip(annotation_pattern, record))
                            for record in record_holder)
        granges = list()
        for record in annotated_holder:
            name = re.search(name_pattern, record['data']).group(1)
            name = name if name else None
            granges.append(GenomicRange(name=name, genome=genome, **record))
        return GenomicRangesList(granges, genome=genome)


def extract_taxid_mapping(mapping, taxid):
    taxid_mapping = dict()
    for key, value in mapping.items():
        try:
            taxid_mapping[key] = value[taxid]
        except KeyError:
            continue
    return taxid_mapping
