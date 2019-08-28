import json
import re
from pathlib import Path
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import namedtuple, defaultdict
from sortedcontainers import SortedKeyList
from .alignment_utils import Alignment


class GenomicException(Exception):
    pass


class ChromosomeNotFoundError(GenomicException):

    def __init__(self, chrom, genome, list_of_chroms):
        self.chrom = chrom
        self.genome = genome
        self.list_of_chroms = list_of_chroms
        super().__init__(f"Chromosome {self.chrom} is not "
                         f"in {self.genome} among these "
                         f"chromosomes: {', '.join(list_of_chroms)}")


class InconsistentChromosomesError(GenomicException):

    def __init__(self, one_chrom, other_chrom):
        self.one_chrom = one_chrom
        self.other_chrom = other_chrom
        super().__init__(f"Inconsistent chromosomes of "
                         f"comparing genomic ranges: "
                         f"{self.one_chrom} and {self.other_chrom}.")


class InconsistentGenomesError(GenomicException):

    def __init__(self, one_genome, other_genome):
        self.one_genome = one_genome
        self.other_genome = other_genome
        super().__init__(f"Inconsistent genomes of "
                         f"comparing genomic ranges: "
                         f"{self.one_genome} and {self.other_genome}.")


class EmptyGenomicRangesListError(GenomicException):

    def __init__(self, range_list_instance):
        self.range_list = range_list_instance
        super().__init__(f"Empty GenomicRangesList instance: {self.range_list}")


class GenomicCoordinatesError(GenomicException):

    def __init__(self, grange, chromosome):
        self.grange = grange
        self.chromosome = chromosome
        super().__init__(f"Coordinates of GenomicRange {self.grange} "
                         f"are out of chromosome with size {self.chromosome.size}")


class LocationNotSpecified(GenomicException):

    def __init__(self, instance):
        self.instance = instance
        super().__init__(f"Sequence file location for {self.instance} "
                         f"is not specified.")


def fasta_reformatter(infile, outfile, total, linewidth=60):
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
    def __init__(self, chrom, start, end, strand='.', name=None,
                 genome=None, sequence_file_loc=None,
                 synteny=None, neighbours=None, **kwargs):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.genome = genome
        self._sequence_file_loc = sequence_file_loc
        self.sequence_header = f"{self.chrom}:{self.start}" \
                               f"-{self.end}{self.strand}"
        self.name = self.sequence_header if name is None else name
        self.connections = {'synteny': (GenomicRangesList([], self.genome)
                                        if synteny is None else synteny),
                            'neighbours': (GenomicRangesList([], self.genome)
                                           if neighbours is None else neighbours)}

    def __repr__(self):
        return f"GenomicRange({self.chrom}, {self.start}, {self.end}, {self.strand})"

    def __str__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}"

    def distance(self, other):
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        if self.genome != other.genome:
            raise InconsistentGenomesError(self.genome, other.genome)
        if self.end < other.start:
            return other.start - self.end
        if other.end < self.start:
            return other.end - self.start
        return 0

    def merge(self, other):
        if self.chrom != other.chrom:
            raise InconsistentChromosomesError(self.chrom, other.chrom)
        if self.genome != other.genome:
            raise InconsistentGenomesError(self.genome, other.genome)
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        return GenomicRange(chrom=self.chrom,
                            start=start,
                            end=end,
                            strand=".",
                            genome=self.genome)

    def align(self, other, **kwargs):
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

    def align_with_connections(self, connection='synteny', **kwargs):
        if connection not in self.connections.keys():
            raise ValueError(f"Connection {connection} not in "
                             f"available list of connections: "
                             f"{self.connections.keys()}.")
        alignment_list = list()
        for grange in self.connections[connection]:
            alignment_list.append(self.align(grange))
        return alignment_list

    def merge_connections(self, connection, distance=0):
        try:
            self.connections[connection] = self.connections[connection].merge(distance)
        except KeyError:
            raise KeyError(f"There is no {connection} connections for the {self}.")

    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(self, fileobj):
        raise NotImplementedError


class AlignedRangePair:

    def __init__(self, query_grange, subject_grange, alignment):
        self.query_grange = query_grange
        self.subject_grange = subject_grange
        self.alignment = alignment

    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(self, fileobj):
        raise NotImplementedError


ChromosomeLocation = namedtuple('ChromosomeLocation', 'size start')


class FastaSeqFile:
    def __init__(self, location):
        self.location = location
        self._chromsizes = dict()
        self._file_obj = None

    @property
    def chromsizes(self):
        if not self._chromsizes:
            if self.location is None:
                raise LocationNotSpecified(self)
            with open(self.location, 'r') as infile:
                chrom, size, start = None, None, None
                line = infile.readline()
                while line:
                    if line.startswith(">"):
                        self._chromsizes[chrom] = ChromosomeLocation(size, start)
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
        if self.location is None:
            raise LocationNotSpecified(self)
        self._file_obj = open(self.location, 'r')

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file_obj.close()

    def get_fasta_by_coord(self, grange, outfile):
        try:
            chrom_size, chrom_start = self.chromsizes[grange.chrom]
        except KeyError:
            raise ChromosomeNotFoundError(grange.chrom,
                                          self.location,
                                          self.chromsizes.keys())
        if grange.end > chrom_size:
            raise GenomicCoordinatesError(grange, self.chromsizes[grange.chrom])
        self._file_obj.seek(chrom_start + grange.start, 0)
        outfile.write(">" + grange.sequence_header + "\n")
        fasta_reformatter(self._file_obj,
                          outfile,
                          total=(grange.end - grange.start),
                          linewidth=60)


class GenomicRangesList(SortedKeyList):
    # TODO: add reading from file (bed, gff, gtf)
    # TODO: parsing orthodb-synteny output
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
                    self[current_self].connections['neighbours'].add(other[current_other])
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

    def connection_mapping(self, mapping, connection):
        for key, value in mapping.items():
            try:
                self.name_mapping[key].connections[connection].append(value)
            except KeyError:
                continue

    @classmethod
    def parse_annotation(cls,
                         fileobj,
                         filetype,
                         genome=None,
                         annotation_pattern=None,
                         name_pattern=None,
                         comment=None,
                         sep=None):
        if filetype not in cls.filetypes:
            raise ValueError(f"filetype {filetype} not one of {cls.filetypes}.")
        if annotation_pattern is None:
            annotation_pattern = cls.annotation_patterns[filetype]
        if name_pattern is None:
            name_pattern = cls.name_patterns[filetype]
        if comment is None:
            comment = cls.comments[filetype]
        if sep is None:
            sep = cls.seps[filetype]
        holder = [dict(zip(annotation_pattern, line.strip().split(sep)))
                  for line in fileobj
                  if line[0] != comment]
        granges = list()
        for record in holder:
            name = re.search(name_pattern, record['data']).group(1)
            name = name if name else None
            granges.append(GenomicRange(name=name, genome=genome, **record))
        return GenomicRangesList(granges, genome=genome)


def translate_orthodb_mapping(self, mapping, grangelist, taxid):
    grange_mapping = defaultdict(list)
    for key, value in mapping.items():
        try:
            map_ids = value[taxid]
            for gene_id in map_ids:
                try:
                    grange_mapping[key].append(grangelist[gene_id])
                except KeyError:
                    continue
        except KeyError:
            continue
    return grange_mapping