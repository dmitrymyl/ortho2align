from collections import namedtuple
from itertools import zip_longest
from tqdm import tqdm
from .genomicranges import BaseGenomicRange, BaseGenomicRangesList


class LiftOverException(Exception):
    """Basic exception for LiftOver."""
    pass


class DeletedInNewError(LiftOverException):
    """Occurs when grange intersects no chains."""

    def __init__(self, grange):
        self.grange = grange
        super().__init__(f'Genomic range {self.grange} intersects none provided chains.')


class PartiallyDeletedInNewError(LiftOverException):
    """Occurs when grange insufficiently intersects one chain."""

    def __init__(self, grange, chain):
        self.grange = grange
        self.chain = chain
        super().__init(f'Genomic range {self.grange} insufficiently intersects one chain: {self.chain}.')


class SplitInNewError(LiftOverException):
    """Occurs when grange insufficiently intersects multiple chains."""

    def __init__(self, grange, chains):
        self.grange = grange
        self.chains = chains
        super().__init__(f'Genomic range {self.grange} insufficiently intersects multiple chains: {self.chains}.')


class DuplicatedInNewError(LiftOverException):
    """Occurs when grange sufficiently intersects multiple chains and that is undesired situation."""

    def __init__(self, grange, chains):
        self.grange = grange
        self.chains = chains
        super().__init__(f'Genomic range {self.grange} sufficiently intersects multiple chains: {self.chains} ' \
                         'and that is undesired situation.')


class InsideTheGapError(LiftOverException):
    """Occurs when grange fully lies within a gap inside a chain."""

    def __init__(self, grange, chain):
        self.grange = grange
        self.chain = chain
        super().__init__(f'Genomic range {self.grange} fully lies within a gap inside a chain {self.chain}.')


class Chain:
    """Convenient class for chain representation when parsing."""
    __slots__ = ('chain', 'score', 'qchrom', 'qsize', 'qstrand', 'qstart',
                 'qend', 'schrom', 'ssize', 'sstrand', 'sstart', 'send', 'name')

    def __init__(self, chain, score, qchrom, qsize, qstrand, qstart,
                 qend, schrom, ssize, sstrand, sstart, send, name):
        self.chain = chain
        self.score = score
        self.qchrom = qchrom
        self.qsize = int(qsize)
        self.qstrand = qstrand
        if qstrand == '-':
            self.qstart = self.qsize - int(qend)
            self.qend = self.qsize - int(qstart)
        elif qstrand in ['+', '.']:
            self.qstart = int(qstart)
            self.qend = int(qend)
        else:
            raise ValueError('The strand of the lift over chain is not one of "+", "-", ".".')
        self.schrom = schrom
        self.ssize = int(ssize)
        self.sstrand = sstrand
        if sstrand == '-':
            self.sstart = self.ssize - int(send)
            self.send = self.ssize - int(sstart)
        elif sstrand in ['+', '.']:
            self.sstart = int(sstart)
            self.send = int(send)
        else:
            raise ValueError('The strand of the lift over chain is not one of "+", "-", ".".')
        self.name = name


class LiftOverChain(BaseGenomicRange):
    """Class for representing one side of liftOver chain.

    A child of BaseGenomicRange class extends its attributes
    with `chromsize`, `blocks`, `gaps`, `sister` via `__slots__`.
    Provides additional methods for liftOver process. Coordinates
    of the chain as genomic range are in the genomic coordinate
    system (start and end are counted from the beginning of the
    chromosome), but internal processes perform switching to the
    liftOver coordinate system (start and end are counted from the
    beginning of the strand, i.e. for '-' strand start and end are
    counted from the end of the chromosome).

    liftOver methods:
        _get_block_indices
        _lift_coordinates

    Properties:
        other_strand
    """
    __slots__ = ('chromsize', 'blocks', 'gaps', 'sister')

    def __init__(self, chrom, start, end, strand, name, chromsize, blocks, gaps, sister=None):
        """Initializes LiftOverChain instance.

        Args:
            chrom (str, int): name of the chromosome.
            start (int): start of the genomic range.
            end (int): end of the genomic range.
            strand (str): strand of the genomic range,
                one of `"+"`, `"-"`, `"."`.
            name (str): name of the genomic range.
            chromsize (int): size of the chromosome.
            blocks (tuple): a tuple containing sizes of
                alignment blocks in a chain.
            gaps (tuple): a tuple containing sizes of
                gaps from this chain.
            sister (LiftOverChain): sister liftOver chain
                (default: None).

        Returns:
            None
        """
        super().__init__(chrom=chrom,
                         start=start,
                         end=end,
                         strand=strand,
                         name=name)
        self.chromsize = chromsize
        self.blocks = blocks
        self.gaps = gaps
        self.sister = sister

    @property
    def init_args(self):
        """Returns __init__ arguments except for sister."""
        return self.chrom, self.start, self.end, self.strand, \
               self.name, self.chromsize, self.blocks, self.gaps

    def __repr__(self):
        """Returns code representation of the liftOver chain."""
        return f"LiftOverChain({', '.join([repr(i) for i in self.init_args])})"

    def __str__(self):
        """Returns tab-separated table representation of the liftOver chain."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t" \
               f"{self.strand}\t{self.chromsize}\t{self.blocks}\t{self.gaps}"

    @property
    def other_strand(self):
        """Returns the opposite strand of the chain's one in a special way.

        The opposite of '+' is '-' and vice versa, but
        the opposite of '.' is '.'

        Raises:
            ValueError: an error in case strand of the chain is not one
            of `'+'`, `'-'`, `'.'`.
        """
        if self.strand == '+':
            return '-'
        elif self.strand == '-':
            return '+'
        elif self.strand == '.':
            return '.'
        else:
            raise ValueError('The strand of the lift over chain is not one of "+", "-", ".".')

    def _get_block_indices(self, grange):
        """Finds indices of start and end blocks that genomic range intersects.

        A private method to be used in the liftOver process. This method
        finds leftmost and rightmost alignments blocks that genomic range
        intersects. Inside performs switching between genomic and liftOver
        coordinate system.

        Args:
            grange (BasicGenomicRange): genomic range to use.

        Returns:
            int, int, int, int: index of the start block, index of
            the end block, offset of genomic range start from the start
            of the start block and offset of genomic range end from the
            start of the end block.

        Raises:
            ValueError: an error in case provided genomic range lies
            outside the chain.
            ValueError: an error in case the strand of the chain is
            not one of `"+"`, `"-"`, `"."`.
            InsideTheGapError: an error in case genomic range lies
            fully inside the gap.
        """
        if grange.end < self.start or grange.start > self.end:
            raise ValueError('Provided grange lies outside this chain.')

        if self.strand == '-':
            current_block_start = self.chromsize - self.end
            grange_start = self.chromsize - grange.end
            grange_end = self.chromsize - grange.start
        elif self.strand in ['+', '.']:
            current_block_start = self.start
            grange_start = grange.start
            grange_end = grange.end
        else:
            raise ValueError('The strand of the lift over chain is not one of "+", "-", ".".') 

        search_for_start = True
        for i, (block_size, gap_size) in enumerate(zip_longest(self.blocks, self.gaps)):
            if search_for_start:  # searcing for the grange_start
                if current_block_start + block_size >= grange_start:  # grange_start inside the block or in the previous gap
                    start_index = i
                    start_offset = grange_start - current_block_start
                    search_for_start = False
                elif gap_size is None: # in practice, should never happen
                    raise ValueError('Provided grange lies outside this chain.')

            if not search_for_start:  # searching for the grange_end
                if (gap_size is None or
                    current_block_start + block_size >= grange_end or 
                    current_block_start + block_size + gap_size >= grange_end):  # grange_end inside the block, the previous gap or downstream the end of the chain
                    end_index = i
                    end_offset = grange_end - current_block_start
                    break
            current_block_start += block_size + gap_size

        if grange_end < current_block_start:  # grange is fully inside one gap
            raise InsideTheGapError(grange, self)
        return start_index, end_index, start_offset, end_offset

    def _lift_coordinates(self, start_index, end_index, start_offset, end_offset):
        """Returns lifted coordinates of the genomic range.

        A private method to be used in the liftOver method. Taken
        the output of the `self.sister._get_block_indices` composes
        lifted coordinates of the genomic range in this chain as follows:
        finds starts of start block and end block by corresponding indices
        and adds corresponding offsets for start and end coordinates.
        Performs switching between genomic and liftOver coordinates.

        Args:
            start_index (int): index of the start block.
            end_index (int): index of the end block.
            start_offset (int): offset of the genomic range from
                the start of the start block.
            end_offset (int): offset of the genomic range from
                the end of the end block.

        Returns:
            int, int: start and end lifted coordinates of the genomic range.

        Raises:
            ValueError: an error in case the strand of the chain is not one
            of "+", "-", ".".
        """
        start_in_chain = sum(self.blocks[:start_index]) + sum(self.gaps[:start_index]) + start_offset
        end_in_chain = sum(self.blocks[:end_index]) + sum(self.gaps[:end_index]) + end_offset
        if self.strand == '-':
            grange_start = self.end - end_in_chain
            grange_end = self.end - start_in_chain
        elif self.strand in ['+', '.']:
            grange_start = self.start + start_in_chain
            grange_end = self.start + end_in_chain
        else:
            raise ValueError('The strand of the lift over chain is not one of "+", "-", ".".')
        return grange_start, grange_end


LiftOverResults = namedtuple('LiftOverResults',
                             'lifted deleted_in_new partially_deleted split_in_new duplicated_in_new inside_the_gap')


class LiftOverChains:
    """Represents liftOver chain file contents.

    Contains two `BaseGenomicRangesList` instances:
    `query_chains` and `subject_chains`, which
    are collections of `LiftOverChain` instances
    representing query and subject sides of chains.
    Instance must be created via `LiftOverChains.parse_chain_file()`.

    Attributes:
        query_chains
        subject_chains

    Basic methods:
        __init__
        __str__
        __repr__

    liftOver methods:
        parse_chain_file
        lift_grange

    Class methods:
        parse_chain_file
    """

    def __init__(self, query_chains, subject_chains):
        """Initializes LiftOverChains instance.

        Args:
            query_chains (BaseGenomicRangesList): query sides of chains.
                A collection of `LiftOverChain` instances.
            subject_chains (BaseGenomicRangesList): subject sides of chains.
                A collection of `LiftOverChain` instances.

        Returns:
            None
        """
        self.query_chains = query_chains
        self.subject_chains = subject_chains

    def __repr__(self):
        """Returns code representation of LiftOverChains instance.

        Relations between sisters are not included.
        """
        return f"LiftOverChains({repr(self.query_chains)}, {repr(self.subject_chains)})"

    def __str__(self):
        """Returns tab-separated table representation of LiftOverChains instance.

        Relations between sisters are not included.
        """
        return f"liftOver chains\nQuery side:\n{str(self.query_chains)}\n" \
               f"Subject side:\n{str(self.subject_chains)}"

    @classmethod
    def parse_chain_file(cls, fileobj, verbose=False):
        """Parses provided liftOver chain file.

        Args:
            fileobj (stream): stream with chain file.
            verbose (bool): whether to show progress
                of parsing via `tqdm` progress bar or
                not (default: False).

        Returns:
            LiftOverChains: parsed liftOver chain file.

        Raises:
            ValueError: in case strand of any chain is not
            one of "+", "-", ".".
            ValueError: in case there are no liftOver chains
            parsed.
        """
        query_chains = list()
        subject_chains = list()
        record_lengths = (3, 1)
        for line in tqdm(fileobj, disable=not verbose):
            record = tuple(line.strip().split())
            if line.startswith('chain'):
                chain_ = Chain(*record)
                blocks = list()
            elif len(record) in record_lengths:
                blocks.append(record)
            elif len(record) == 0:
                sizes = tuple(int(item[0]) for item in blocks)
                dqs = tuple(int(item[1]) for item in blocks[:-1])
                dss = tuple(int(item[2]) for item in blocks[:-1])

                query_chain = LiftOverChain(chrom=chain_.qchrom,
                                            start=chain_.qstart,
                                            end=chain_.qend,
                                            strand=chain_.qstrand,
                                            name=chain_.name,
                                            blocks=sizes,
                                            gaps=dqs,
                                            chromsize=chain_.qsize)
                subject_chain = LiftOverChain(chrom=chain_.schrom,
                                              start=chain_.sstart,
                                              end=chain_.send,
                                              strand=chain_.sstrand,
                                              name=chain_.name,
                                              blocks=sizes,
                                              gaps=dss,
                                              chromsize=chain_.ssize)

                query_chain.sister = subject_chain
                subject_chain.sister = query_chain
                query_chains.append(query_chain)
                subject_chains.append(subject_chain)
            else:
                if len(query_chains) == 0 or len(subject_chains) == 0:
                    raise ValueError('Provided file is not in the correct format of chain file.')
        return cls(query_chains=BaseGenomicRangesList(query_chains),
                   subject_chains=BaseGenomicRangesList(subject_chains))

    def lift_grange(self, grange, min_ratio, origin='query', allow_duplications=False):
        """Lifts one genomic range.

        Args:
            grange (BaseGenomicRange): genomic range to be lifted.
            min_ratio (int, float): minimal fraction of the genomic
                range overlapping with a chain to consider the latter
                sufficient for lifting.
            origin (str): which genome the genomic range
                came from. Must be one of "query", "subject"
                (default: "query").
            allow_duplications (bool): if True, will allow
                sufficient overlap with more than one chain
                (default: False).

        Returns:
            BaseGenomicRangesList: a list of lifted genomic ranges.

        Raises:
            ValueError: in case `origin` is not one of "query", "subject".
            DeletedInNewError: in case genomic range intersects no chains.
            PartiallyDeletedInNewError: in case genomic range
            insufficiently intersects one chain.
            SplitInNewError: in case genomic range insufficiently
            intersects more than one chain.
            DuplicatedInNewError: in case genomic range sufficiently
            intersects more than one chain and `allow_duplicatons`
            set to `False`.
            InsideTheGapError: in case genomic range falls within the
            gap of all sufficiently intersecting chains.
        """
        if origin == 'query':
            query_chains = self.query_chains
        elif origin == 'subject':
            query_chains = self.subject_chains
        else:
            raise ValueError('origin is not one of "query", "subject".')

        intersecting_chains = grange.find_neighbours(query_chains, distance=0)

        if len(intersecting_chains) == 0:
            raise DeletedInNewError(grange)

        sufficient_chains = BaseGenomicRangesList([chain_
                                                   for chain_ in intersecting_chains
                                                   if grange.calc_fraction(chain_) >= min_ratio])

        if len(intersecting_chains) == 1 and len(sufficient_chains) == 0:
            raise PartiallyDeletedInNewError(grange, intersecting_chains[0])
        if len(intersecting_chains) > 1 and len(sufficient_chains) == 0:
            raise SplitInNewError(grange, intersecting_chains)
        if len(sufficient_chains) > 1 and not allow_duplications:
            raise DuplicatedInNewError(grange, sufficient_chains)

        lifted_granges = list()
        for i, query_chain in enumerate(sufficient_chains):
            try:
                lifted_start, lifted_end = query_chain.sister \
                                                      ._lift_coordinates(*query_chain._get_block_indices(grange))
            except ValueError as e:
                raise e
            except InsideTheGapError as e:
                if len(sufficient_chains) - i == 1 and len(lifted_granges) == 0:
                    raise e
                continue

            if query_chain.strand == grange.strand:
                lifted_strand = query_chain.sister.strand
            else:
                lifted_strand = query_chain.sister.other_strand
            lifted_chrom = query_chain.sister.chrom
            if grange._name is not None:
                lifted_name = grange.name
            else:
                lifted_name = None
            used_args = {'chrom', 'start', 'end', 'strand', 'name'}
            kwargs = {attr: getattr(grange, attr)
                      for attr in set(grange.init_args) - used_args}
            lifted_granges.append(grange.__class__(chrom=lifted_chrom,
                                                   start=lifted_start,
                                                   end=lifted_end,
                                                   strand=lifted_strand,
                                                   name=lifted_name,
                                                   **kwargs))
        return BaseGenomicRangesList(lifted_granges)

    def lift_granges(self, granges_list, min_ratio, origin='query', allow_duplications=False):
        """Lifts genomic ranges list.

        Args:
            granges_list (BaseGenomicRangeList or children): genomic range
                to be lifted.
            min_ratio (int, float): minimal fraction of the genomic
                range overlapping with a chain to consider the latter
                sufficient for lifting. The suggested value for intraspecies
                liftOver is 0.95, for interspecies is 0.1.
            origin (str, default "query"): which genome the genomic range
                came from. Must be one of "query", "subject".
            allow_duplications (bool, default False): if True, will allow
                sufficient overlap with more than one chain.

        Returns:
            tuple of 6 granges_list.__class__ instances: A list of lifted
            genomic ranges, a list of genomic ranges deleted in new genome,
            a list of genomic ranges partially deleted in new genome,
            a list of genomic ranges split in new genome,
            a list of genomic ranges duplicated in new genome
            and a list of genomic ranges found inside the gap
            of liftOver chains.

        Raises:
            ValueError: An error in case *origin* is not one of *"query"*,
                *"subject"*.
        """

        if origin not in ['query', 'subject']:
            raise ValueError('origin is not one of "query", "subject".')
        lifted = []
        deleted_in_new = []
        partially_deleted = []
        split_in_new = []
        duplicated_in_new = []
        inside_the_gap = []

        for grange in granges_list:
            try:
                lifted_granges = self.lift_grange(grange,
                                                  min_ratio,
                                                  origin,
                                                  allow_duplications)
                lifted += [new_grange for new_grange in lifted_granges]
            except DeletedInNewError:
                deleted_in_new.append(grange)
            except PartiallyDeletedInNewError:
                partially_deleted.append(grange)
            except SplitInNewError:
                split_in_new.append(grange)
            except DuplicatedInNewError:
                duplicated_in_new.append(grange)
            except InsideTheGapError:
                inside_the_gap.append(grange)
        used_args = {'collection', }
        kwargs = {attr: getattr(granges_list, attr)
                  for attr in set(granges_list.init_args) - used_args}
        return LiftOverResults(*map(lambda new_list: granges_list.__class__(new_list, **kwargs),
                                    (lifted,
                                     deleted_in_new,
                                     partially_deleted,
                                     split_in_new,
                                     duplicated_in_new,
                                     inside_the_gap)))
