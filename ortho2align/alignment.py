import matplotlib.pyplot as plt


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
        (bool) True or False, result of comparison
        of a and b with side function.
    Raises:
        ValueError in case side is not one of 'g',
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
        (str): "+" if both strands are concordant, "-" otherwise.
            "." is equal to "+".
    """
    _alignment_strand = alignment_strand in ['+', '.']
    if hsp_strand == _alignment_strand:
        return "+"
    else:
        return "-"


class HSP:
    """
    High-scoring pair from BLAST alignment representing
    an alignment block.
    """
    orientation_dict = {(True, True): 'direct',
                        (True, False): 'reverse'}

    def __init__(self, qstart, qend, sstart, send, **kwargs):
        """Initializes HSP class.

        Args:
            qstart (int): start of HSP in query sequence.
            qend (int): end of HSP in query sequence.
            sstart (int): start of HSP in subject sequence.
            send (int): end of HSP in subject sequence.
            score (int, float): alignment score of HSP.
            qchrom (int, str): chromosome of query sequence
                (default: None).
            schrom (int, str): chromosome of subject sequence
                (default: None).
        """
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.kwargs = kwargs
        self.qstrand = self.qend - self.qstart > 0
        self.sstrand = self.send - self.sstart > 0
        self.orientation = self.orientation_dict.get((self.qstrand,
                                                      self.sstrand))

    def __repr__(self):
        kwargs_repr = self.kwargs.__repr__()
        params_list = [str(i)
                       for i in (self.qstart,
                                 self.qend,
                                 self.sstart,
                                 self.send,
                                 self.score)] + ["**" + kwargs_repr]
        return f"HSP({', '.join(params_list)})"

    def __str__(self):
        return f"q {self.qstart}:{self.qend} s {self.sstart}:{self.send} " \
               f"score {self.score} {self.orientation}"

    def __eq__(self, other):
        """Equality magic method."""
        return (self.qstart == other.qstart and
                self.qend == other.qend and
                self.sstart == other.sstart and
                self.send == other.send and
                self.score == other.score)

    def precede(self, other):
        """Infers HSP precedence.

        Infers whether self HSP precedes other HSP
        in both query and subject sequence.

        Returns:
            True if self precedes other, False otherwise.

        Raises:
            ValueError in case of inconsistent orientations.
            Precedence is defined only for HSPs of same
            orientation (except for None orientation).
        """
        if (self.orientation != other.orientation and
            (self.orientation is not None) and
            (other.orientation is not None)):
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
            self: HSP instance;
            other: HSP instance;
            gapopen (int): penalty for gap opening;
            gapextend (int): penalty for gap extension.

        Returns:
            An int/float or Inf in case none HSP precedes
            the other one.

        Raises:
            ValueError in case of inconsistent orientations.
            Distance is defined only for HSPs of same
            orientation (except for None orientation).
        """
        if (self.orientation != other.orientation and
            (self.orientation is not None) and
            (other.orientation is not None)):
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
        """Returns a dict representation."""
        return self.__dict__

    @classmethod
    def from_dict(cls, dict_):
        """Recovers from a dict representation."""
        keys = ['qstart', 'qend', 'sstart', 'send', 'score']
        return cls(*[dict_.get(key, None) for key in keys],
                   **dict_.get('kwargs', dict()))

    def copy(self):
        """Returns a copy of an instance.

        All dynamic data in HSPVertex will be lost!
        """
        return HSP(self.qstart,
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
            None if inplace, modified copy of hsp otherwise.
        """
        if not inplace:
            hsp = hsp.copy()

        if type_ == 'qleft':
            fraction = (hsp.qend - point) / (hsp.qend - hsp.qstart)
            hsp.sstart = round((point - hsp.qstart) *
                               (hsp.send - hsp.sstart) /
                               (hsp.qend - hsp.qstart) +
                               hsp.sstart)
            hsp.qstart = point
        elif type_ == 'qright':
            fraction = (point - hsp.qstart) / (hsp.qend - hsp.qstart)
            hsp.send = round((point - hsp.qstart) *
                             (hsp.send - hsp.sstart) /
                             (hsp.qend - hsp.qstart) +
                             hsp.sstart)
            hsp.qend = point
        elif type_ == 'sleft':
            qpoint = round((point - hsp.sstart) *
                           (hsp.qend - hsp.qstart) /
                           (hsp.send - hsp.sstart) +
                           hsp.qstart)
            if hsp.orientation == 'direct':
                fraction = (hsp.send - point) / (hsp.send - hsp.sstart)
                hsp.qstart = qpoint
                hsp.sstart = point
            else:
                fraction = (hsp.sstart - point) / (hsp.sstart - hsp.send)
                hsp.qend = qpoint
                hsp.send = point
        elif type_ == 'sright':
            qpoint = round((point - hsp.sstart) *
                           (hsp.qend - hsp.qstart) /
                           (hsp.send - hsp.sstart) +
                           hsp.qstart)
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
    """An extension of HSP class for handling graph
    operations (contains a list of neighbours).
    """

    def __init__(self, qstart, qend, sstart, send, score, **kwargs):
        """Inits HSPVertex class.

        Args:
            qstart (int): start of HSP in query sequence.
            qend (int): end of HSP in query sequence.
            sstart (int): start of HSP in subject sequence.
            send (int): end of HSP in subject sequence.
            qchrom (int, str): chromosome of query sequence
                (default: None).
            schrom (int, str): chromosome of subject sequence
                (default: None).
        """
        super().__init__(qstart, qend, sstart, send, score, **kwargs)
        self.next_vertices = []
        self.total_score = float("Inf")
        self.best_prev = None

    def __repr__(self):
        kwargs_repr = self.kwargs.__repr__()
        params_list = [str(_)
                       for _ in (self.qstart,
                                 self.qend,
                                 self.sstart,
                                 self.send,
                                 self.score)] + ["**" + kwargs_repr]
        return f"HSPVertex({', '.join(params_list)})"

    def __str__(self):
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


def is_float(string):
    """Checks whether a string represent a float number."""
    try:
        _ = float(string)
        return True
    except ValueError:
        return False


def numberize(string):
    """String to number if possible.

    Transforms string to float or to
    integer if possible.

    Returns:
        int, float or str instance.
    """
    if is_float(string):
        a = float(string)
        if a.is_integer():
            return int(a)
        else:
            return a
    else:
        return string


class Alignment:
    """BLAST alignment of two sequences.

    Represents a BLAST alignment of two sequences
   (NOT an alignment of a sequence against a database!).
    Contains a list of HSPs and other additional information.
    Allows to read BLAST alignment from file produced with
    `-outfmt 7`, plot an alignment map, read from and write to
    dictionary representation (hence, JSON) and find best
    transcripts based on some strategies.

    Attributes:
        _all_HSPs (list): list of all HSPs related
            to the alignment.
        HSPs (list): list of HSP instances that were
            filtered by some criterion.
        qlen (int): length of query sequence.
        slen (int): length of subject sequence.
        replace_dict (dict): field names to replace
            for convenience.

    """
    replace_dict = {'q. start': 'qstart',
                    'q. end': 'qend',
                    's. start': 'sstart',
                    's. end': 'send',
                    'query length': 'qlen',
                    'subject length': 'slen'}

    def __init__(self,
                 HSPs,
                 qlen=None,
                 slen=None,
                 filtered_HSPs=None):
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
        self.HSPs = HSPs.copy() if filtered_HSPs is None else filtered_HSPs
        self.qlen = (max(self._all_HSPs,
                         key=lambda hsp: hsp.qend).qend + 1
                     if qlen is None else qlen)
        self.slen = (max(max(self._all_HSPs,
                             key=lambda hsp: hsp.sstart).sstart + 1,
                         max(self._all_HSPs,
                             key=lambda hsp: hsp.send).send + 1)
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

    def plot_alignment(self, qleft=None, qright=None, sleft=None, sright=None):
        """Plots alignment map."""
        for hsp in self.HSPs:
            plt.plot([hsp.qstart, hsp.qend],
                     [hsp.sstart, hsp.send],
                     color='black')
        for hsp in set(self._all_HSPs) - set(self.HSPs):
            plt.plot([hsp.qstart, hsp.qend],
                     [hsp.sstart, hsp.send],
                     color='grey')
        plt.xlim(qleft, qright)
        plt.ylim(sleft, sright)
        plt.xlabel("Query")
        plt.ylabel("Subject")

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
            An Alignment instance.

        Raises:
            ValueError in case field 'score'
            is not found in alignment fields.
        """
        line = file_object.readline()
        if not line.startswith('#'):
            raise ValueError("Provided file_object is not in accepted file format.")

        while not line.startswith("# Fields:"):
            line = file_object.readline()

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
                   hsp.kwargs.get('qlen'),
                   hsp.kwargs.get('slen'),
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

    def reset_filter_by_score(self):
        """Resets HSP filtering by score."""
        self.HSPs = [hsp for hsp in self._all_HSPs]
        self.filtered = False

    def cut_coordinates(self, qleft=None, qright=None, sleft=None, sright=None):
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
            New Alignment instance with HSPs within given
            coordinates.
        """
        hsps = [hsp.copy() for hsp in self._all_HSPs]
        if qleft is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.qstart < qleft <= hsp.qend:
                    HSP._process_boundary(hsp, qleft, 'qleft')
                    new_hsps.append(hsp)
                    pass
                elif qleft <= hsp.qstart:
                    new_hsps.append(hsp)
            hsps = new_hsps
        if qright is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.qstart <= qright < hsp.qend:
                    HSP._process_boundary(hsp, qright, 'qright')
                    new_hsps.append(hsp)
                elif hsp.qend <= qright:
                    new_hsps.append(hsp)
            hsps = new_hsps
        if sleft is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.sstart < sleft <= hsp.send:
                    HSP._process_boundary(hsp, sleft, 'sleft')
                    new_hsps.append(hsp)
                    pass
                elif sleft <= hsp.sstart and hsp.orientation == 'direct':
                    new_hsps.append(hsp)
                elif hsp.send < sleft <= hsp.sstart:
                    HSP._process_boundary(hsp, sleft, 'sleft')
                    new_hsps.append(hsp)
                    pass
                elif sleft <= hsp.send and hsp.orientation == 'reverse':
                    new_hsps.append(hsp)
            hsps = new_hsps
        if sright is not None:
            new_hsps = list()
            for hsp in hsps:
                if hsp.sstart <= sright < hsp.send:
                    HSP._process_boundary(hsp, sright, 'sright')
                    new_hsps.append(hsp)
                    pass
                elif hsp.send <= sright and hsp.orientation == 'direct':
                    new_hsps.append(hsp)
                elif hsp.send <= sright < hsp.sstart:
                    HSP._process_boundary(hsp, sright, 'sright')
                    new_hsps.append(hsp)
                    pass
                elif hsp.sstart <= sright and hsp.orientation == 'reverse':
                    new_hsps.append(hsp)
            hsps = new_hsps
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
            Alignment instance with new HSP coordinates
            and new qlen and slen.
        """
        # TODO: deal with reverse-derived sequences and their representation.
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
            a new list of HSPVertex instances.
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

    def _find_all_transcripts(self, vertices, orientation):
        """
        Build all possible transcripts from given list
        of HSPs so that HSPs follow each other on both
        query and subject sequences and do not overlap
        with each other.

        Args:
            vertices (list): a list of HSP instances
                sorted by query start and query end.

        Returns:
            list of all possible transcripts.
        """
        if len(vertices) == 0:
            return []
        vertices = self._build_hsp_graph(vertices, orientation)
        stack = list()
        transcripts = list()
        stack.append([vertices[0]])
        while stack:
            if stack[-1][-1].next_vertices:
                nascent_transcript = stack.pop()
                for neighbour, weight in nascent_transcript[-1].next_vertices:
                    stack.append(nascent_transcript + [neighbour])
            else:
                transcripts.append(Transcript(stack.pop()[1:-1], self))
        return transcripts

    def get_all_transcripts(self):
        """Finds all possible transcripts from HSPs.

        Finds all possible transcripts consisting of
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
            dictionary containing lists of transcripts
            in direct and reverse orientation.
        """
        oriented_groups = self._split_orientations()
        return {key: self._find_all_transcripts(group, key)
                for key, group in oriented_groups.items()}

    def _find_best_transcript(self, vertices, orientation):
        """Finds best transcript in given HSP list.

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
            A Transcript instance with the best score.
        """
        if len(vertices) == 0:
            return Transcript([], self)
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
        return Transcript(best_hsps[1:][::-1], self, -score)

    def get_best_transcripts(self):
        """Finds best transcripts.

        Finds best transcripts in direct
        and reverse HSP orientation with
        dynamic programming approach (i.e.
        minimizing gap penalties between HSPs.)

        Returns:
            dict with best direct and reverse
            Transcript instance under 'direct'
            and 'reverse' keys, respectively.
        """
        oriented_groups = self._split_orientations()
        return {key: self._find_best_transcript(group, key)
                for key, group in oriented_groups.items()}

    def best_transcript(self):
        """Finds the best transcript based on its score.

        First, finds best transcripts in direct and
        reverse orientations and then choses the one
        with the biggest score.

        Returns:
            Transcript: transcript of the biggest score.
        """
        transcripts = self.get_best_transcripts()
        if transcripts['direct'].score > transcripts['reverse'].score:
            return transcripts['direct']
        else:
            return transcripts['reverse']


class Transcript:
    """Transcript combined from non-overlapping
    sequential HSPs.

    Attributes:
        HSPs (list): a list of HSPVertex instances
            representing a sequence of homologous
            regions in query and subject sequences;
        alignment (Alignment): corresponding
            `Alignment` instance;
        score (int, float): a transcript score
            computed as HSPs scores minus gap
            penalties (see `set_score` and
            `HSP.distance` methods for definition).
    """

    def __init__(self, HSPs, alignment, score=None):
        """Inits Transcript class.

        Args:
            HSPs (list): a list of HSPs representing the
                transcript.
            alignment (GenomicRangesAlignment): corresponding `Alignment`
                instance.
        """
        self.HSPs = HSPs
        self.alignment = alignment
        self._score = score

    @property
    def score(self):
        """Calculates transcript score.

        The transcript score is based on alignment
        score and gap penalties defined by distance
        function.

        Returns:
            Transcript score.
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
        return f"Transcript({repr(self.HSPs)})"

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

    def to_bed12(self, mode='list'):
        """
        Turns transcript into BED12 representation of both
        query side and subject side.

        Args:
            mode (str): how to return representation. If 'list' then
            list of lists, if 'str' then complete BED12 record.

        Returns:
            (list of two lists or two strs): BED12 representation.
        """
        # TODO: test
        if len(self.HSPs) == 0:
            raise ValueError('No HSPs were found in the Transcript.')
        Qchrom = self.alignment.qchrom
        QchromStart = min(min(self.HSPs, key=lambda i: i.qstart),
                          min(self.HSPs, key=lambda i: i.qend))
        QchromEnd = max(max(self.HSPs, key=lambda i: i.qstart),
                        max(self.HSPs, key=lambda i: i.qend))
        Qname = self.alignment.qname
        Qscore = 1000
        Qstrand = nxor_strands(self.HSPs[0].qstrand, self.alignment.qstrand)
        QthickStart = QchromStart
        QthickEnd = QchromEnd
        QitemRgb = 0
        QblockCount = len(self.HSPs)
        QblockSizes = [abs(hsp.qend - hsp.qstart) for hsp in self.HSPs]
        if Qstrand == "+":
            QblockStarts = [hsp.qstart - QchromStart for hsp in self.HSPs]
        else:
            QblockStarts = [hsp.qend - QchromStart for hsp in self.HSPs]

        Schrom = self.alignment.schrom
        SchromStart = min(min(self.HSPs, key=lambda i: i.sstart),
                          min(self.HSPs, key=lambda i: i.send))
        SchromEnd = max(max(self.HSPs, key=lambda i: i.sstart),
                        max(self.HSPs, key=lambda i: i.send))
        Sname = self.alignment.Sname
        Sscore = 1000
        Sstrand = nxor_strands(self.HSPs[0].sstrand, self.alignment.sstrand)
        SthickStart = SchromStart
        SthickEnd = SchromEnd
        SitemRgb = 0
        SblockCount = len(self.HSPs)
        SblockSizes = [abs(hsp.send - hsp.sstart) for hsp in self.HSPs]
        if Sstrand == "+":
            SblockStarts = [hsp.sstart - SchromStart for hsp in self.HSPs]
        else:
            SblockStarts = [hsp.Send - SchromStart for hsp in self.HSPs]

        q_side = [Qchrom, QchromStart, QchromEnd, Qname, Qscore,
                  Qstrand, QthickStart, QthickEnd, QitemRgb, QblockCount,
                  QblockSizes, QblockStarts]
        s_side = [Schrom, SchromStart, SchromEnd, Sname, Sscore,
                  Sstrand, SthickStart, SthickEnd, SitemRgb, SblockCount,
                  SblockSizes, SblockStarts]
        if mode == 'list':
            return [q_side, s_side]
        elif mode == 'str':
            return ["\t".join([','.join(item) if isinstance(item, list) else item for item in q_side]),
                    "\t".join([','.join(item) if isinstance(item, list) else item for item in s_side])]
        else:
            raise ValueError('mode not one of ["list", "str"].')

    def plot_transcript(self, color='red', link_color='blue'):
        """Plots transcript with corresponding alignment.

        First plots corresponding alignment in black, then
        transcripts's HSPs with solid lines and connections
        between them with dashed lines.

        Args:
            color (str): a color name to plot transcript's
                HSPs;
            link_color (str): a color name to plot connections
                between transcript's HSPs.

        Returns:
            None
        """
        self.alignment.plot_alignment()
        for hsp in self.HSPs:
            plt.plot([hsp.qstart, hsp.qend],
                     [hsp.sstart, hsp.send],
                     color=color)
        for i in range(len(self.HSPs) - 1):
            plt.plot([self.HSPs[i].qend, self.HSPs[i + 1].qstart],
                     [self.HSPs[i].send, self.HSPs[i + 1].sstart],
                     color=link_color,
                     linestyle='dashed')
