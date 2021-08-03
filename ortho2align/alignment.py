from .utils import numberize


def compare(a, b, side='g'):
    """Calculates boolean function of two values.

    Args:
        a: the first argument to compare
            with the second argument.
        b: the second argument to compare
            with the first argument.
        side (str): boolean function code which
            is used to compare a with b. The accepted
            values are:
            'g': greater, a > b
            'l': less, a < b
            'eq': equal, a == b
            'neq': not equal, a != b
            'geq': greater or equal, a >= b
            'leq': less or equal, a <= b
    Returns:
        bool: True or False, result of comparison
            of a and b with side function.
    Raises:
        ValueError: in case side is not one of 'g',
            'l', 'eq', 'neq', 'geq', 'leq'.
    """
    if side == 'g':
        return a > b
    elif side == 'l':
        return a < b
    elif side == 'eq':
        return a == b
    elif side == 'neq':
        return a != b
    elif side == 'geq':
        return a >= b
    elif side == 'leq':
        return a <= b
    else:
        raise ValueError('Incorrect side argument.')


def nxor_strands(hsp_strand, alignment_strand):
    """Returns NXOR of hsp strand and alignment strand.

    Args:
        hsp_strand (boolean): HSP.qstrand or HSP.sstrand.
        qlignment_strand (str): Alignment.qstrand or Alignment.sstrand.
            One of "+", "-", ".".

    Returns:
        str: "+" if both strands are concordant, "-" otherwise.
            "." is equal to "+".
    """
    _alignment_strand = alignment_strand in ['+', '.']
    if hsp_strand == _alignment_strand:
        return "+"
    else:
        return "-"


class HSP:
    """High-scoring pair from BLAST alignment representing an alignment block.

    Attributes:
        qstart (int): start of HSP in query sequence.
        qend (int): end of HSP in query sequence.
        sstart (int): start of HSP in subject sequence.
        send (int): end of HSP in subject sequence.
        score (int, float): alignment score of HSP.
        kwargs (dict): kwargs from ``init`` method.
        qstrand (bool): True if query strand is "+", False otherwise.
        sstrand (bool): True if subject strand is "+", False otherwise.
        orientation (str): 'direct' if both strands are "+" and 'reverse' otherwise.

    Class attributes:
        orientation_dict (dict): holds information for HSP orientation.
    """
    orientation_dict = {(True, True): 'direct',
                        (True, False): 'reverse'}
    __slots__ = ('qstart', 'qend', 'sstart', 'send', 'score',
                 'kwargs', 'qstrand', 'sstrand', 'orientation')

    def __init__(self, qstart, qend, sstart, send, score, **kwargs):
        """Initializes HSP class.

        Args:
            qstart (int): start of HSP in query sequence.
            qend (int): end of HSP in query sequence.
            sstart (int): start of HSP in subject sequence.
            send (int): end of HSP in subject sequence.
            score (int, float): alignment score of HSP.
            kwargs (dict): any other keyword arguments.
        """
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.score = score
        self.kwargs = kwargs
        self.qstrand = self.qend - self.qstart > 0
        self.sstrand = self.send - self.sstart > 0
        self.orientation = self.orientation_dict.get((self.qstrand,
                                                      self.sstrand))

    def __repr__(self):
        """Representation of HSP instance."""
        kwargs_repr = self.kwargs.__repr__()
        params_list = [str(i)
                       for i in (self.qstart,
                                 self.qend,
                                 self.sstart,
                                 self.send,
                                 self.score)] + ["**" + kwargs_repr]
        return f"HSP({', '.join(params_list)})"

    def __str__(self):
        """String representation of HSP instance."""
        return f"q {self.qstart}:{self.qend} s {self.sstart}:{self.send} " \
               f"score {self.score} {self.orientation}"

    def __eq__(self, other):
        """Equality magic method."""
        return (self.qstart == other.qstart
                and self.qend == other.qend
                and self.sstart == other.sstart
                and self.send == other.send
                and self.score == other.score)

    def precede(self, other):
        """Infers HSP precedence.

        Infers whether self HSP precedes other HSP
        in both query and subject sequence.

        Args:
            other (HSP): another instance of HSP class to
                infer precedence.

        Returns:
            bool: True if self precedes other, False otherwise.

        Raises:
            ValueError: in case of inconsistent orientations.
                Precedence is defined only for HSPs of same
                orientation (except for None orientation).
        """
        if self.orientation != other.orientation \
           and (self.orientation is not None) \
           and (other.orientation is not None):
            raise ValueError(f"HSPs have inconsistent orientations: "
                             f"{self.orientation} for {self} and "
                             f"{other.orientation} for {other}.")
        qprecede = self.qend < other.qstart
        if self.orientation == 'direct' or other.orientation == 'direct':
            sprecede = self.send < other.sstart
        else:
            sprecede = self.send > other.sstart
        if qprecede and sprecede:
            return True
        return False

    def distance(self, other, gapopen=5, gapextend=2):
        """Calculates distance between two HSPs.

        Distance between two HSPs of same orientation
        is calculated as sum of distances in query and
        subject sequence multiplied by gap extension
        penalty added with gap opening penalty and
        subtracted with the score of the following HSP.
        In case none HSP precedes the other one (overlap
        or different, non-sequential placement) the
        distance is set to infinity.

        Args:
            other: HSP instance.
            gapopen (int): penalty for gap opening.
            gapextend (int): penalty for gap extension.

        Returns:
            int, float: An int/float distance or Inf in case none HSP precedes
                the other one.

        Raises:
            ValueError: in case of inconsistent orientations.
                Distance is defined only for HSPs of same
                orientation (except for None orientation).
        """
        if self.orientation != other.orientation \
           and (self.orientation is not None) \
           and (other.orientation is not None):
            raise ValueError(f"HSPs have inconsistent orientations: "
                             f"{self.orientation} for {self} and "
                             f"{other.orientation} for {other}.")
        if self.precede(other):
            qdist = other.qstart - self.qend
            if self.orientation == 'direct' or other.orientation == 'direct':
                sdist = other.sstart - self.send
            else:
                sdist = self.send - other.sstart
            return gapopen + gapextend * (qdist + sdist) - other.score
        elif other.precede(self):
            qdist = self.qstart - other.qend
            if self.orientation == 'direct' or other.orientation == 'direct':
                sdist = self.sstart - other.send
            else:
                sdist = other.send - self.sstart
            return gapopen + gapextend * (qdist + sdist) - self.score
        else:  # doesn't distinguish between overlap and different placement
            return float("Inf")

    def to_dict(self):
        """Returns a dict representation.

        Returns:
            dict: a dict representation of the instance.
        """
        return {'qstart': self.qstart,
                'qend': self.qend,
                'sstart': self.sstart,
                'send': self.send,
                'score': self.score,
                'kwargs': self.kwargs}

    @classmethod
    def from_dict(cls, dict_):
        """Recovers from a dict representation.

        Dict must contain the following keys:
        'qstart', 'qend', 'sstart', 'send', 'score'.

        Args:
            dict_ (dict): a dictionary representation of
                HSP instance generated with `HSP.to_dict` method.

        Returns:
            HSP: HSP instance.
        """
        keys = ['qstart', 'qend', 'sstart', 'send', 'score']
        return cls(*[dict_.get(key, None) for key in keys],
                   **dict_.get('kwargs', dict()))

    def copy(self):
        """Returns a shallow copy of the instance.

        All dynamic data in HSPVertex will be lost!

        Returns:
            HSP: a shallow copy of the instance.
        """
        return self.__class__(self.qstart,
                              self.qend,
                              self.sstart,
                              self.send,
                              self.score,
                              **self.kwargs)

    @staticmethod
    def _process_boundary(hsp, point, type_, inplace=True):
        """Cuts HSP to a certain boundary.

        In case an HSP intersects a boundary, it is cut
        and its score is reduced proportionally to the
        share of the cut.

        Args:
            hsp (HSP): an HSP instance to modify;
            point (int, float): a coordinate to be set
                new boundary;
            type_ (str): a type of point. One of 'qleft',
                'qright', 'sleft', 'sright';
            inplace (bool): whether to modify HSP in-place
                or return a modified copy (default: True)

        Returns:
            HSP, None: None if inplace, modified copy of hsp otherwise.
        """
        if not inplace:
            hsp = hsp.copy()

        if type_ == 'qleft':
            fraction = (hsp.qend - point) / (hsp.qend - hsp.qstart)
            hsp.sstart = round((point - hsp.qstart)
                               * (hsp.send - hsp.sstart)
                               / (hsp.qend - hsp.qstart)
                               + hsp.sstart)
            hsp.qstart = point
        elif type_ == 'qright':
            fraction = (point - hsp.qstart) / (hsp.qend - hsp.qstart)
            hsp.send = round((point - hsp.qstart)
                             * (hsp.send - hsp.sstart)
                             / (hsp.qend - hsp.qstart)
                             + hsp.sstart)
            hsp.qend = point
        elif type_ == 'sleft':
            qpoint = round((point - hsp.sstart)
                           * (hsp.qend - hsp.qstart)
                           / (hsp.send - hsp.sstart)
                           + hsp.qstart)
            if hsp.orientation == 'direct':
                fraction = (hsp.send - point) / (hsp.send - hsp.sstart)
                hsp.qstart = qpoint
                hsp.sstart = point
            else:
                fraction = (hsp.sstart - point) / (hsp.sstart - hsp.send)
                hsp.qend = qpoint
                hsp.send = point
        elif type_ == 'sright':
            qpoint = round((point - hsp.sstart)
                           * (hsp.qend - hsp.qstart)
                           / (hsp.send - hsp.sstart)
                           + hsp.qstart)
            if hsp.orientation == 'direct':
                fraction = (point - hsp.sstart) / (hsp.send - hsp.sstart)
                hsp.qend = qpoint
                hsp.send = point
            else:
                fraction = (point - hsp.send) / (hsp.sstart - hsp.send)
                hsp.qstart = qpoint
                hsp.sstart = point
        hsp.score *= fraction

        if not inplace:
            return hsp


class HSPVertex(HSP):
    """An extension of HSP class for handling graph operations.

    Contains a list of next vertices needed for dynamic programming.

    Attributes:
        qstart (int): start of HSP in query sequence.
        qend (int): end of HSP in query sequence.
        sstart (int): start of HSP in subject sequence.
        send (int): end of HSP in subject sequence.
        score (int, float): alignment score of HSP.
        kwargs (dict): any other keyword arguments.
        next_vertices (list): a list of the following vertices.
        total_score (float): a total score of chain ending on this HSP.
        best_prev (HSPVertex): the best preceding HSP in terms of score.
    """
    __slots__ = ('next_vertices', 'total_score', 'best_prev')

    def __init__(self, qstart, qend, sstart, send, score, **kwargs):
        """Inits HSPVertex class.

        Args:
            qstart (int): start of HSP in query sequence.
            qend (int): end of HSP in query sequence.
            sstart (int): start of HSP in subject sequence.
            send (int): end of HSP in subject sequence.
            score (int, float): alignment score of HSP.
            kwargs (dict): any other keyword arguments.
        """
        super().__init__(qstart, qend, sstart, send, score, **kwargs)
        self.next_vertices = []
        self.total_score = float("Inf")
        self.best_prev = None

    def __repr__(self):
        """Representation of ``HSPVertex`` instance."""
        kwargs_repr = self.kwargs.__repr__()
        params_list = [str(_)
                       for _ in (self.qstart,
                                 self.qend,
                                 self.sstart,
                                 self.send,
                                 self.score)] + ["**" + kwargs_repr]
        return f"HSPVertex({', '.join(params_list)})"

    def __str__(self):
        """String epresentation of ``HSPVertex`` instance."""
        next_vertices_str = "\t" + "\n\t".join(" ".join(str(i).split(" ")[:-1]) for i in self.next_vertices)
        return f"q {self.qstart}:{self.qend} s {self.sstart}:{self.send} score " \
               f"{self.score} {self.orientation}\n" + next_vertices_str

    def _relax(self, other, weight):
        """A relax method for dynamic programming.

        If the current path from the begin vertex
        to self vertex is more expensive than the path
        that goes through the other vertex, a new cheaper
        path from the begin to self is set.

        Args:
            self (HSPVertex): a vertex to be relaxed,
                i.e. to minimize score of path from
                begin vertex to the vertex.
            other (HSPVertex): a vertex preceding
                and connected to self.
            weight (int, float): weight of the edge
                connecting other and self vertices.

        Returns:
            None
        """
        new_score = other.total_score + weight
        if self.total_score > new_score:
            self.total_score = new_score
            self.best_prev = other


class Alignment:
    """BLAST alignment of two sequences.

    Represents a BLAST alignment of two sequences
   (NOT an alignment of a sequence against a database!).
    Contains a list of HSPs and other additional information.
    Allows to read BLAST alignment from file produced with
    `-outfmt 7`, plot an alignment map, read from and write to
    dictionary representation (hence, JSON) and find best
    alignment chains based on some strategies.

    Attributes:
        _all_HSPs (list): list of all HSPs related
            to the alignment.
        HSPs (list): list of HSP instances that were
            filtered by some criterion.
        qlen (int): length of query sequence.
        slen (int): length of subject sequence.
        replace_dict (dict): field names to replace
            for convenience.

    Class attributes:
        replace_dict (dict): keys are BLAST+ table headers,
            values are corresponding ``HSP`` attributes.

    """
    replace_dict = {'q. start': 'qstart',
                    'q. end': 'qend',
                    's. start': 'sstart',
                    's. end': 'send',
                    'query length': 'qlen',
                    'subject length': 'slen'}
    __slots__ = ('_all_HSPs', 'HSPs', 'qlen', 'slen', 'filtered')

    def __init__(self,
                 HSPs,
                 qlen=None,
                 slen=None,
                 filtered_HSPs=None,
                 filtered=None):
        """Inits Alignment class.

        Args:
            HSPs (list): list of HSP instances;
            qlen (int): query sequence length
                (default: None);
            slen (int): subject sequence length
                (default: None);
            qchrom (str, int): query chromosome name (default: None).
            schrom (str, int): subject chromosome name (default: None).
            qname (str): query sequence name (default: None).
            sname (str): subject sequence name (default: None).
            qstrand (str): query sequence name.
                One of "+", "-", "." (default: ".").
            sstrand (str): subject sequence name.
                One of "+", "-", "." (default: ".").
            filtered_HSPs (list): list of HSP instances
                that were filtered by some criterion
               (default: None).

        """
        self._all_HSPs = HSPs
        self.HSPs = [hsp.copy() for hsp in HSPs] if filtered_HSPs is None else filtered_HSPs
        self.qlen = (max(self._all_HSPs,
                         key=lambda hsp: hsp.qend,
                         default=HSPVertex(0, 0, 0, 0, 0)).qend + 1
                     if qlen is None else qlen)
        self.slen = (max(max(self._all_HSPs,
                             key=lambda hsp: hsp.sstart,
                             default=HSPVertex(0, 0, 0, 0, 0)).sstart + 1,
                         max(self._all_HSPs,
                             key=lambda hsp: hsp.send,
                             default=HSPVertex(0, 0, 0, 0, 0)).send + 1)
                     if slen is None else slen)
        self.filtered = False if filtered_HSPs is None else True

    def __repr__(self):
        return f"Alignment({self.HSPs.__repr__()}, {self.qlen}, {self.slen})"

    def __str__(self):
        return f"Alignment of {self.qlen} and {self.slen}\n" + \
               "\n".join(hsp.__str__() for hsp in self.HSPs)

    def __eq__(self, other):
        """Equality magic method."""
        if len(self._all_HSPs) != len(other._all_HSPs):
            return False
        if self.qlen != other.qlen or self.slen != other.slen:
            return False
        return all([i == k for i, k in zip(self._all_HSPs, other._all_HSPs)])

    @property
    def chain_class(self):
        """Corresponding alignment chain class.

        Returns:
            type: `AlignmentChain`.
        """
        return AlignmentChain

    @classmethod
    def from_file_blast(cls,
                        file_object,
                        start_type=1,
                        end_type='inclusive',
                        **kwargs):
        """File parser.

        Takes file object that corresponds to
        a BLAST alignment of two sequences
        produced with `-outfmt 7` key and
        parses it into an Alignment object.

        Args:
            file_object: file object.
            qchrom (str, int): query chromosome name (default: None).
            schrom (str, int): subject chromosome name (default: None).
            qname (str): query sequence name (default: None).
            sname (str): subject sequence name (default: None).
            qstrand (str): query sequence name.
                One of "+", "-", "." (default: ".").
            sstrand (str): subject sequence name.
                One of "+", "-", "." (default: ".").
            start_type (int): chromosome start coordinate.
                One of 0, 1 (default: 1).
            end_type (str): chromosome end inclusion type.
                One of 'inclusive', 'exclusive' (default: 'inclusive').

        Returns:
            Alignment: An `Alignment` instance.

        Raises:
            ValueError: in case field 'score'
                is not found in alignment fields.
        """
        line = file_object.readline()
        if not line.startswith('#'):
            raise ValueError("Provided file_object is not in accepted file format.")

        for _ in range(100):
            line = file_object.readline()
            if line.startswith("# Fields:"):
                break
        else:
            return cls([],
                       **kwargs)

        fields = line.lstrip("# Fields: ").rstrip("\n").split(", ")
        if 'score' not in fields:
            raise ValueError("'score' is not found in alignment fields. "
                             "Please provide an alignment file with score "
                             "field (raw score, not bit score).")
        fields = [Alignment.replace_dict.get(item, item)
                  for item in fields] + ['qchrom', 'schrom']
        hit_number = int(file_object.readline().split(" ")[1])
        HSPs = list()
        for i in range(hit_number):
            hsp = file_object.readline().strip().split("\t")
            hsp = [numberize(item) for item in hsp]
            hsp = HSPVertex(**dict(zip(fields, hsp)))
            if start_type == 1:
                hsp.qstart -= 1
                hsp.qend -= 1
                hsp.sstart -= 1
                hsp.send -= 1
            if end_type == 'inclusive':
                hsp.qend += 1
                hsp.send += 1
            HSPs.append(hsp)
        return cls(HSPs,
                   qlen=hsp.kwargs.get('qlen'),
                   slen=hsp.kwargs.get('slen'),
                   **kwargs)

    def to_dict(self):
        """Transforms instance to dict representation."""
        keys = ['HSPs',
                'qlen',
                'slen',
                'filtered_HSPs',
                'filtered']
        HSPs = [hsp.to_dict() for hsp in self._all_HSPs]
        filtered_HSPs = [hsp.to_dict() for hsp in self.HSPs]
        return dict(zip(keys, [HSPs,
                               self.qlen,
                               self.slen,
                               filtered_HSPs,
                               self.filtered]))

    @classmethod
    def from_dict(cls, dict_):
        """Restores object from dict representation."""
        HSPs = [HSPVertex.from_dict(item)
                for item in dict_.get('HSPs', [])]
        qlen = dict_.get('qlen')
        slen = dict_.get('slen')
        filtered = dict_.get('filtered')
        if dict_.get('filtered_HSPs') is None:
            filtered_HSPs = None
        else:
            filtered_HSPs = [HSPVertex.from_dict(item)
                             for item in dict_.get('filtered_HSPs')]
        return cls(HSPs, qlen, slen, filtered_HSPs, filtered)

    def filter_by_score(self, score, side='g'):
        """Filters HSPs by score.

        All HSPs being stored in `self.HSPs`,
        the method puts those HSPs in `self.HSPs`
        that have score compared to `score` in a
        manner specified by `side` argument. As default (side='g'),
        HSP.score have to be greater than `score`.

        Args:
            self: Alignment instance;
            score (int): score to filter by.
            side (str): a side to compare score and HSP
                scores by. One of:
                'g' for greater (HSP.score > score)
                'l' for less (HSP.score < score)
                'eq' for equal (HSP.score == score)
                'neq' for not equal (HSP.score != score)
                'geq' for greater or equal (HSP.score >= score)
                'leq' for less or equal (HSP.score <= score)
                Default: 'g'.

        Returns:
            None
        """
        self.HSPs = [hsp
                     for hsp in self.HSPs
                     if compare(hsp.score, score, side)]
        self.filtered = True

    def filter_by_function_score(self, function, score, side='g'):
        """Filter HSPs by function value of score.

        All HSPs being stored in `self.HSPs`,
        the method puts those HSPs in `self.HSPs`
        that have function value of score compared to `score` in a
        manner specified by `side` argument. As default (side='g'),
        function(HSP.score) have to be greater than `score`.

        Args:
            self: Alignment instance;
            function (callable): a function to apply to
                HSP.score.
            score (int): score to filter by.
            side (str): a side to compare score and function results
                by. One of:
                'g' for greater (HSP.score > score)
                'l' for less (HSP.score < score)
                'eq' for equal (HSP.score == score)
                'neq' for not equal (HSP.score != score)
                'geq' for greater or equal (HSP.score >= score)
                'leq' for less or equal (HSP.score <= score)
                Default: 'g'.

        Returns:
            None
        """
        self.HSPs = [hsp
                     for hsp in self.HSPs
                     if compare(function(hsp.score), score, side)]
        self.filtered = True

    def filter_by_array(self, array, score, side='g'):
        """Filters HSPs by corresponding values in array.

        There is a corresponding value in `array` for
        every HSP in `self.HSPs`. The method filters only
        those HSPs, which corresponding value is greater/less/e.t.c.
        than the `score` (the type of comparison is specified by `side`).

        Args:
            self: Alignment instance;
            array (list, tuple, e.t.c.): a sequence of corresponding
                values to `self.HSPs`.
            score (int): score to filter by.
            side (str): a side to compare score and function results
                by. One of:
                'g' for greater (HSP.score > score)
                'l' for less (HSP.score < score)
                'eq' for equal (HSP.score == score)
                'neq' for not equal (HSP.score != score)
                'geq' for greater or equal (HSP.score >= score)
                'leq' for less or equal (HSP.score <= score)
                Default: 'g'.
        Raises:
            ValueError: in case `len(array) != len(self.HSPs)`.
        """
        if len(array) != len(self.HSPs):
            raise ValueError(f'Length of provided array is not equal to the '
                             f'amount of HSPs: {len(array)} vs {len(self.HSPs)}.')
        self.HSPs = [hsp
                     for hsp, value in zip(self.HSPs, array)
                     if compare(value, score, side)]

    def filter_by_bool_array(self, array):
        """Filters HSPs by True/False values in the array.

        Args:
            self: Alignment instance;
            array (collection of bools): a sequence of corresponding
                bool values to `self.HSPs`.
        Raises:
            ValueError: in case `len(array) != len(self.HSPs)`.
        """
        if len(array) != len(self.HSPs):
            raise ValueError(f'Length of provided array is not equal to the '
                             f'amount of HSPs: {len(array)} vs {len(self.HSPs)}.')
        self.HSPs = [hsp
                     for hsp, indicator in zip(self.HSPs, array)
                     if indicator]

    def reset_filter_by_score(self):
        """Resets HSP filtering by score."""
        self.HSPs = [hsp.copy() for hsp in self._all_HSPs]
        self.filtered = False

    def _cut_hsps(self, qleft=None, qright=None, sleft=None, sright=None):
        """Leaves only HSPs in certain range of coordinates.

        Given query and subject sequence boundaries only
        those HSPs are kept that lay within those boundaries.
        If an HSP intersects a boundary, it is cut and its
        score is reduced proportionally to the share of the cut.
        If a boundary is left None than it is equal to +Inf
        or -Inf for left and right boundaries, respectively.

        Args:
            qleft (int, float): query sequence left boundary
                (default: None);
            qright (int, float): query sequence right boundary
                (default: None);
            sleft (int, float): subject sequence left boundary
                (default: None);
            sright (int, float): subjecct sequence right boundary
                (default: None);

        Returns:
            list: list of cut hsps.
        """
        hsps = [hsp.copy() for hsp in self._all_HSPs]
        if qleft is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.qstart < qleft <= hsp.qend:
                    HSPVertex._process_boundary(hsp, qleft, 'qleft')
                    new_hsps.append(hsp)
                elif qleft <= hsp.qstart:
                    new_hsps.append(hsp)
            hsps = new_hsps
        if qright is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.qstart <= qright < hsp.qend:
                    HSPVertex._process_boundary(hsp, qright, 'qright')
                    new_hsps.append(hsp)
                elif hsp.qend <= qright:
                    new_hsps.append(hsp)
            hsps = new_hsps
        if sleft is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.sstart < sleft <= hsp.send:
                    HSPVertex._process_boundary(hsp, sleft, 'sleft')
                    new_hsps.append(hsp)
                elif sleft <= hsp.sstart and hsp.orientation == 'direct':
                    new_hsps.append(hsp)
                elif hsp.send < sleft <= hsp.sstart:
                    HSPVertex._process_boundary(hsp, sleft, 'sleft')
                    new_hsps.append(hsp)
                elif sleft <= hsp.send and hsp.orientation == 'reverse':
                    new_hsps.append(hsp)
            hsps = new_hsps
        if sright is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.sstart <= sright < hsp.send:
                    HSPVertex._process_boundary(hsp, sright, 'sright')
                    new_hsps.append(hsp)
                elif hsp.send <= sright and hsp.orientation == 'direct':
                    new_hsps.append(hsp)
                elif hsp.send <= sright < hsp.sstart:
                    HSPVertex._process_boundary(hsp, sright, 'sright')
                    new_hsps.append(hsp)
                elif hsp.sstart <= sright and hsp.orientation == 'reverse':
                    new_hsps.append(hsp)
        return hsps

    def cut_coordinates(self, qleft=None, qright=None, sleft=None, sright=None):
        """Creates a new ``Alignment`` instance with HSPs restricted by given coordinates.

        Given query and subject sequence boundaries only
        those HSPs are kept that lay within those boundaries.
        If an HSP intersects a boundary, it is cut and its
        score is reduced proportionally to the share of the cut.
        If a boundary is left None than it is equal to +Inf
        or -Inf for left and right boundaries, respectively.

        Args:
            qleft (int, float): query sequence left boundary
                (default: None);
            qright (int, float): query sequence right boundary
                (default: None);
            sleft (int, float): subject sequence left boundary
                (default: None);
            sright (int, float): subjecct sequence right boundary
                (default: None);

        Returns:
            Alignment: a new alignment instance with cut HSPs.
        """
        hsps = self._cut_hsps(qleft=qleft, qright=qright, sleft=sleft, sright=sright)
        return Alignment(hsps, self.qlen, self.slen)

    def set_zero(self, qzero, szero):
        """Sets new query and subject sequence beginnings.

        Given coordinates of the current query and subject
        sequence beginnnings in the new coordinates,
        transforms all HSP coordinates accordingly. Returns
        new Alignment instance, filtration and graph data
        will be lost. Use before filtering and graph building.

        Args:
            qzero (int): the location of the current query
                sequence beginning in the new coordinates.
            szero (int): the location of the current subject
                sequence beginning in the new coordinates.

        Returns:
            Alignment: `Alignment` instance with new HSP coordinates
            and new qlen and slen.
        """
        # Not dealing with reversed sequences, but ok.
        qlen = self.qlen + qzero
        slen = self.slen + szero
        HSPs = list()
        for hsp in self._all_HSPs:
            new_hsp = hsp.copy()
            new_hsp.qstart += qzero
            new_hsp.qend += qzero
            new_hsp.sstart += szero
            new_hsp.send += szero
            HSPs.append(new_hsp)
        return Alignment(HSPs, qlen, slen)

    def _split_orientations(self):
        """Splits HSPs in two groups by orientation."""
        sorted_hsps = sorted(self.HSPs,
                             key=lambda hsp: (hsp.orientation,
                                              hsp.qstart,
                                              hsp.sstart))
        oriented_groups = {'direct': [],
                           'reverse': []}
        for hsp in sorted_hsps:
            oriented_groups[hsp.orientation].append(hsp)
        return oriented_groups

    def _build_hsp_graph(self, vertices, orientation):
        """Builds graph from HSP list.

        Taken sorted list of HSPVertex instances
        adds succeeding HSPs in `next_vertices` attribute so
        as these HSPs are the closest to the given
       (only one if it does not overlap any other HSP,
        more than one otherwise). Check the code and
        `HSP.precede` for complete definition.

        Args:
            vertices (list): a list of HSPVertex instances
                of one direction sorted by qstart.

        Returns:
            list: a new list of HSPVertex instances.
        """
        # Add begin and end nodes.
        if orientation == 'direct':
            begin = HSPVertex(*[0 for _ in range(5)])
            end = HSPVertex(self.qlen,
                            self.qlen,
                            self.slen,
                            self.slen,
                            0)
        else:
            begin = HSPVertex(0, 0, self.slen, self.slen, 0)
            end = HSPVertex(self.qlen, self.qlen, 0, 0, 0)
        vertices = [begin] + vertices + [end]
        # Join all those vertex pairs where one
        # vertex precedes the other one.
        for i in range(len(vertices)):
            vertices[i].next_vertices = list()
            vertices[i].best_prev = None
            vertices[i].total_score = 0 if i == 0 else float("Inf")
            for j in range(i + 1, len(vertices)):
                if vertices[i].precede(vertices[j]):
                    vertices[i].next_vertices.append([vertices[j],
                                                      vertices[i].distance(vertices[j])])
        # Leave only those edges that connect each vertex
        # with its direct preceder(s).
        for vertex in vertices:
            i = 0
            k = 1
            while k < len(vertex.next_vertices):
                if vertex.next_vertices[i][0].precede(vertex.next_vertices[k][0]):
                    vertex.next_vertices = vertex.next_vertices[:k]
                else:
                    k += 1

        return vertices

    def _find_all_chains(self, vertices, orientation):
        """
        Build all possible alignment chains from given list
        of HSPs so that HSPs follow each other on both
        query and subject sequences and do not overlap
        with each other.

        Args:
            vertices (list): a list of HSP instances
                sorted by query start and query end.

        Returns:
            list: list of all possible alignment chains.
        """
        if len(vertices) == 0:
            return []
        vertices = self._build_hsp_graph(vertices, orientation)
        stack = list()
        chains = list()
        stack.append([vertices[0]])
        while stack:
            if stack[-1][-1].next_vertices:
                nascent_chain = stack.pop()
                for neighbour, weight in nascent_chain[-1].next_vertices:
                    stack.append(nascent_chain + [neighbour])
            else:
                chains.append(self.chain_class(stack.pop()[1:-1], self))
        return chains

    def get_all_chains(self):
        """Finds all possible alignment chains from HSPs.

        Finds all possible alignment chains consisting of
        consecutive non-overlapping HSPs in both
        query and subject sequences in two orientations:
        * direct, i.e. HSP is aimed from start towards
            end in both query and subject sequences;
        * reverse, i.e. from start to end in query
            sequence but from end to start
            in subject sequence.

        Note:
            Very greedy -- exponential on number of HSPs.
            Use only with few HSPs.

        Returns:
            dictionary containing lists of alignment chains
            in direct and reverse orientation.
        """
        oriented_groups = self._split_orientations()
        return {key: self._find_all_chains(group, key)
                for key, group in oriented_groups.items()}

    def _find_best_chain(self, vertices, orientation):
        """Finds best alignment chain in given HSP list.

        Utilizes dynamic programming algorithm of finding
        the shortest path between two graph vertices to
        get a sequence of HSPs with lowest total distance score.
        Distance score is defined by `HSP.distance` method when
        building the graph with `_build_hsp_graph` method.

        Args:
            vertices (list): list of HSPVertex instances
                produced by `_split_orientations` method
               (i.e. of one orientation and sorted by qstart
                and sstart);
            orientation (str): one of 'direct', 'reverse'.

        Returns:
            An AlignmentChain instance with the best score.
        """
        if len(vertices) == 0:
            return self.chain_class([], self)
        vertices = self._build_hsp_graph(vertices, orientation)
        for vertex in vertices:
            for neighbour, weight in vertex.next_vertices:
                neighbour._relax(vertex, weight)
        current_vertex = vertices[-1]
        score = current_vertex.total_score
        best_hsps = list()
        while current_vertex.best_prev is not None:
            best_hsps.append(current_vertex)
            current_vertex = current_vertex.best_prev
        return self.chain_class(best_hsps[1:][::-1], self, -score)

    def get_best_chains(self):
        """Finds best alignment chains.

        Finds best alignment chains in direct
        and reverse HSP orientation with
        dynamic programming approach (i.e.
        minimizing gap penalties between HSPs.)

        Returns:
            dict: dict with best direct and reverse
                ``AlignmentChain`` instance under 'direct'
                and 'reverse' keys, respectively.
        """
        oriented_groups = self._split_orientations()
        return {key: self._find_best_chain(group, key)
                for key, group in oriented_groups.items()}

    def best_chain(self):
        """Finds the best alignment chain based on its score.

        First, finds best alignment chains in direct and
        reverse orientations and then choses the one
        with the biggest score.

        Returns:
            AlignmentChain: alignment chain of the biggest score.
        """
        chains = self.get_best_chains()
        if chains['direct'].score > chains['reverse'].score:
            return chains['direct']
        else:
            return chains['reverse']


class AlignmentChain:
    """AlignmentChain combined from non-overlapping
    sequential HSPs.

    Attributes:
        HSPs (list): a list of HSPVertex instances
            representing a sequence of homologous
            regions in query and subject sequences;
        alignment (Alignment): corresponding
            `Alignment` instance;
        score (int, float): a alignment chain score
            computed as HSPs scores minus gap
            penalties (see `set_score` and
            `HSP.distance` methods for definition).
    """
    __slots__ = ('HSPs', 'alignment', '_score')

    def __init__(self, HSPs, alignment, score=None):
        """Inits AlignmentChain class.

        Args:
            HSPs (list): a list of HSPs representing the
                alignment chain.
            alignment (GenomicRangesAlignment): corresponding `Alignment`
                instance.
        """
        self.HSPs = HSPs
        self.alignment = alignment
        self._score = score

    @property
    def score(self):
        """Calculates alignment chain score.

        The alignment chain score is based on alignment
        score and gap penalties defined by distance
        function.

        Returns:
            float: AlignmentChain score.
        """
        if self._score is None:
            if self.HSPs:
                score = self.HSPs[0].score
                for i in range(len(self.HSPs) - 1):
                    score -= self.HSPs[i].distance(self.HSPs[i + 1])
                self._score = score
            else:
                self._score = -float('inf')
        return self._score

    def __str__(self):
        header = "qstart\tqend\tsstart\tsend"
        lines = ['\t'.join(str(_)
                           for _ in [hsp.qstart,
                                     hsp.qend,
                                     hsp.sstart,
                                     hsp.send])
                 for hsp in self.HSPs]
        return "\n".join([header] + lines)

    def __repr__(self):
        return f"AlignmentChain({repr(self.HSPs)})"

    def __eq__(self, other):
        if len(self.HSPs) != len(other.HSPs):
            return False
        return all([i == k
                    for i, k in zip([self.HSPs + [self.alignment,
                                                  self.score]],
                                    [other.HSPs + [other.alignment,
                                                   other.score]])])

    def to_dict(self):
        keys = ['HSPs', 'alignment', 'score']
        HSPs = [hsp.to_dict() for hsp in self.HSPs]
        alignment = self.alignment.to_dict()
        return dict(zip(keys, [HSPs,
                               alignment,
                               self.score]))

    @classmethod
    def from_dict(cls, dict_):
        HSPs = [HSPVertex.from_dict(item)
                for item in dict_.get('HSPs', [])]
        alignment = Alignment.from_dict(dict_.get('alignment'))
        score = dict_.get('score')
        return cls(HSPs, alignment, score)
