import re
from tqdm.auto import tqdm
from .genomicranges import GenomicRange, GenomicRangesList


class ParserException(Exception):
    """Base parser exception class."""
    pass


class IncorrectStrand(ParserException):
    """Exception for incorrect strand value."""
    def __init__(self, strand):
        """Initializes an IncorrectStrand instance.

        Args:
            strand (str): an incorrect strand value.
        """
        self.strand = strand

    def __str__(self):
        return f"Incorrect strand value: {self.strand}."


class IncorrectItemRgb(ParserException):
    """Exception for incorrect itemRgb value in a bed12 file."""
    def __init__(self, itemRgb):
        """Initializes an IncorrectItemRgb instance.

        Args:
            itemRgb (str): an incorrect itemRgb value.
        """
        self.itemRgb = itemRgb

    def __str__(self):
        return f"Incorrenct itemRgb value: {self.itemRgb}"


class IncorrectBlockSizes(ParserException):
    """Exception for incorrect blockSizes value in a bed12 file."""
    def __init__(self, blockSizes):
        """Initializes an IncorrectBlockSizes instance.

        Args:
            blockSizes (str): an incorrect blockSizes value.
        """
        self.blockSizes = blockSizes

    def __str__(self):
        return f"Incorrect blockSizes value: {self.blockSizes}"


class InconsistentBlockSizes(IncorrectBlockSizes):
    """Exception for inconsistent blockSizes value in a bed12 file.

    The exception is raised in case the amount of block sizes is not
    equal to the blockCount value.
    """
    def __init__(self, blockSizes, blockCount):
        """Initializes an InconsistentBlockSizes instance.

        Args:
            blockSizes (str): blockSizes value.
            blockCount (str): blockCount value.
        """
        self.blockSizes = blockSizes
        self.blockCount = blockCount

    def __str__(self):
        return f"blockSizes field '{self.blockSizes}' is inconsistent with blockCount field '{self.blockCount}'."


class IncorrectBlockStarts(ParserException):
    """Exception for incorrect blockStarts value in a bed12 file."""
    def __init__(self, blockStarts):
        """Initializes an IncorrectBlockStarts instance.
        
        Args:
            blockStarts (str): blockStarts value.
        """
        self.blockStarts = blockStarts

    def __str__(self):
        return f"Incorrect blockStarts value: {self.blockStarts}"


class InconsistentBlockStarts(IncorrectBlockStarts):
    """Exception for inconsistent blockStarts value in a bed12 file.

    The exception is raised in case the amount of block starts is not
    equal to the blockCount value.
    """
    def __init__(self, blockStarts, blockCount):
        """Initializes an InconsistentBlockStarts instance.

        Args:
            blockStarts (str): blockStarts value.
            blockCount (str): blockCount value.
        """
        self.blockStarts = blockStarts
        self.blockCount = blockCount

    def __str__(self):
        return f"blockStarts field '{self.blockStarts}' is inconsistent with blockCount field '{self.blockCount}'."


class IncorrectPhase(ParserException):
    """Exception for incorrect phase value in a gtf/gff file."""
    def __init__(self, phase):
        """Initializes an IncorrectPhase instance.

        Args:
            phase (str): phase value.s
        """
        self.phase = phase

    def __str__(self):
        return f"Incorrect phase value: {self.phase}"


class IncorrectGTFAttrs(ParserException):
    """Exception for incorrect attributes value in a gtf file."""
    def __init__(self, attrs):
        """Initializes an IncorrectGTFAttrs instance.

        Args:
            attrs (str): attributes value.
        """
        self.attrs = attrs

    def __str__(self):
        return f"Incorrect GTF attributes: {self.attrs}"


class IncorrectGFFAttrs(ParserException):
    """Exception for incorrect attributes value in a gtf file."""
    def __init__(self, attrs):
        """Initializes an IncorrectGFFAttrs instance.

        Args:
            attrs (str): attributes value.
        """
        self.attrs = attrs

    def __str__(self):
        return f"Incorrect GFF attributes: {self.attrs}"


class InconsistentBED3(ParserException):
    """Exception for an inconsistent line in a bed3 file.

    Raised when a number of fields is not equal to 3 in a
    particular line.
    """
    def __init__(self, line):
        """Initializes InconsistentBED3 instance.

        Args:
            line (str): inconsistent line.
        """
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED3 file format in the following line: {self.line}"


class InconsistentBED6(ParserException):
    """Exception for an inconsistent line in a bed6 file.

    Raised when a number of fields is not equal to 5 in a
    particular line.
    """
    def __init__(self, line):
        """Initializes InconsistentBED6 instance.

        Args:
            line (str): inconsistent line.
        """
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED6 file format in the following line: {self.line}"


class InconsistentBED12(ParserException):
    """Exception for an inconsistent line in a bed12 file.

    Raised when a number of fields is not equal to 12 in a
    particular line.
    """
    def __init__(self, line):
        """Initializes InconsistentBED12 instance.

        Args:
            line (str): inconsistent line.
        """
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED12 file format in the following line: {self.line}"


class InconsistentGTF(ParserException):
    """Exception for an inconsistent line in a gtf file.

    Raised when a number of fields is not equal to 9 in a
    particular line.
    """
    def __init__(self, line):
        """Initializes InconsistentGTF instance.

        Args:
            line (str): inconsistent line.
        """
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with GTF file format in the following line: {self.line}"


class InconsistentGFF(ParserException):
    """Exception for an inconsistent line in a GFF file.

    Raised when a number of fields is not equal to 9 in a
    particular line.
    """
    def __init__(self, line):
        """Initializes InconsistentGFF instance.

        Args:
            line (str): inconsistent line.
        """
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with GFF file format in the following line: {self.line}"


class IncorrectLine(ParserException):
    """An exception for a line that cannot be parsed for selected file format."""
    def __init__(self, line, line_no):
        """Initializes an IncorrectLine instance.

        Args:
            line (str): an incorrect line.
            line_no (int): a number of incorrect line (starts with 1).
        """
        self.line = line
        self.line_no = line_no

    def __str__(self):
        return f"Line No. {self.line_no} is incorrect for the selected file format: {self.line}"


class EmptyAnnotation(ParserException):
    """An exception for a file that empty or consists of only the comments."""
    def __str__(self):
        return "Provided annotation file is empty or consist only of comment lines."


class UnrecognizedFormat(ParserException):
    """An exception for an unrecognized file format."""
    def __init__(self, fields):
        self.fields = fields

    def __str__(self):
        return f"Unrecognized annotation format with {self.fields} fields."


def parse_bed_start(field):
    """Parses bed start filed.

    Args:
        field (str): bed start field.

    Returns:
        int: integer start value.
    """
    return int(field)


def parse_bed_end(field):
    """Parses bed end filed.

    Args:
        field (str): bed end field.

    Returns:
        int: integer end value.
    """
    return int(field)


def parse_score(field):
    """Parses bed score filed.

    In case score is an integer, returns integer,
    same goes for a float. If score value is not
    a number, returns string.

    Args:
        field (str): bed score field.

    Returns:
        int, float or str: parsed score value.
    """
    try:
        a = float(field)
        if a.is_integer():
            return int(a)
        else:
            return a
    except ValueError:
        return field


def parse_strand(field):
    """Parses strand field.

    Args:
        field (str): strand field.

    Returns:
        str: parsed strand value.

    Raises:
        IncorrectStrand: in case strand is not
            one of '.', '+', '-'.
    """
    if field in ('.', '+', '-'):
        return field
    raise IncorrectStrand(field)


def parse_itemRgb(field):
    """Parses bed itemRgb field..

    Args:
        field (str): itemRgb field.

    Returns:
        str or tuple: str if field is '0' or '.',
            tuple of 3 ints if field is not empty.

    Raises:
        IncorrectItemRgb: in case itemRgb field
            is not a 3-int csv or cannot be parsed.
    """
    try:
        if field in ('0', '.'):
            return field
        record = tuple(int(x) for x in field.split(','))
        if len(record) != 3:
            raise IncorrectItemRgb(field)
        return record
    except Exception:
        raise IncorrectItemRgb(field)


def parse_blockSizes(field):
    """Parses bed blockSizes field.

    Args:
        field (str): blockSizes field.

    Returns:
        tuple: of ints of block sizes.

    Raises:
        IncorrectBlockSizes: in case field
            cannot be parsed.
    """    
    try:
        record = tuple(int(i) for i in field.split(',') if i != '')
        return record
    except Exception:
        raise IncorrectBlockSizes(field)


def parse_blockStarts(field):
    """Parses bed blockStarts field.

    Args:
        field (str): blockStarts field.

    Returns:
        tuple: of ints of block stars.

    Raises:
        IncorrectBlockStars: in case field
            cannot be parsed.
    """  
    try:
        record = tuple(int(i) for i in field.split(',') if i != '')
        return record
    except Exception:
        raise IncorrectBlockStarts(field)


bed_fields = {'chrom': str,
              'start': parse_bed_start,
              'end': parse_bed_end,
              'name': str,
              'score': parse_score,
              'strand': parse_strand,
              'thickStart': parse_bed_start,
              'thickEnd': parse_bed_end,
              'itemRgb': parse_itemRgb,
              'blockCount': int,
              'blockSizes': parse_blockSizes,
              'blockStarts': parse_blockStarts}
bed_comment = '#'


def bed3_parser(fileobj, verbose=False, sequence_file_path=None):
    """Parses bed3 file into GenomicRangesList instance.

    Args:
        fileobj (file): an opend file stream.
        verbose (bool): if True, reports parsing progress (default: False).
        sequence_file_path (str): path to the corresponding genome file (default: None).

    Returns:
        GenomicRangesList: a list of genomic ranges.

    Raises:
        InconsistentBED3: in case a line contains a number of fields
            other than 3.
        IncorrectLine: in case a line cannot be parsed.
    """
    granges = list()
    colnames = tuple(bed_fields.keys())[:3]
    dtypes = tuple(bed_fields.values())[:3]
    comment = bed_comment

    for line_no, line in tqdm(enumerate(fileobj), disable=not verbose):
        if line.startswith(comment) or line.startswith('track'):
            continue
        try:
            data = tuple(line.strip().split('\t'))
            if len(data) != 3:
                raise InconsistentBED3(line)
            record = {key: func(item)
                      for key, func, item in zip(colnames, dtypes, data)}
            granges.append(GenomicRange(**record))
        except Exception:
            raise IncorrectLine(line, line_no + 1)

    return GenomicRangesList(granges, sequence_file_path)


def bed6_parser(fileobj, verbose=False, sequence_file_path=None):
    """Parses bed6 file into GenomicRangesList instance.

    Args:
        fileobj (file): an opend file stream.
        verbose (bool): if True, reports parsing progress (default: False).
        sequence_file_path (str): path to the corresponding genome file (default: None).

    Returns:
        GenomicRangesList: a list of genomic ranges.

    Raises:
        InconsistentBED6: in case a line contains a number of fields
            other than 6.
        IncorrectLine: in case a line cannot be parsed.
    """
    granges = list()
    colnames = tuple(bed_fields.keys())[:6]
    dtypes = tuple(bed_fields.values())[:6]
    comment = bed_comment

    for line_no, line in tqdm(enumerate(fileobj), disable=not verbose):
        if line.startswith(comment) or line.startswith('track'):
            continue
        try:
            data = tuple(line.strip().split('\t'))
            if len(data) != 6:
                raise InconsistentBED6(line)
            record = {key: func(item)
                      for key, func, item in zip(colnames, dtypes, data)}
            granges.append(GenomicRange(**record))
        except Exception:
            raise IncorrectLine(line, line_no + 1)

    return GenomicRangesList(granges, sequence_file_path)


def bed12_parser(fileobj, verbose=False, sequence_file_path=None):
    """Parses bed12 file into GenomicRangesList instance.

    Args:
        fileobj (file): an opend file stream.
        verbose (bool): if True, reports parsing progress (default: False).
        sequence_file_path (str): path to the corresponding genome file (default: None).

    Returns:
        GenomicRangesList: a list of genomic ranges.

    Raises:
        InconsistentBED12: in case a line contains a number of fields
            other than 12.
        InconsistentBlockSizes: in case of inconsistent blockSizes field.
        InconsistentBlockStarts: in case of inconsistent blockStarts field.
        IncorrectLine: in case a line cannot be parsed.
    """
    granges = list()
    colnames = bed_fields.keys()
    dtypes = bed_fields.values()
    comment = bed_comment

    for line_no, line in tqdm(enumerate(fileobj), disable=not verbose):
        if line.startswith(comment) or line.startswith('track'):
            continue
        try:
            data = tuple(line.strip().split('\t'))
            if len(data) != 12:
                raise InconsistentBED12(line)
            record = {key: func(item)
                      for key, func, item in zip(colnames, dtypes, data)}

            if len(record['blockSizes']) != record['blockCount']:
                raise InconsistentBlockSizes(record['blockSizes'], record['blockCount'])
            if len(record['blockStarts']) != record['blockCount']:
                raise InconsistentBlockStarts(record['blockStarts'], record['blockCount'])
            granges.append(GenomicRange(**record))
        except Exception:
            raise IncorrectLine(line, line_no + 1)

    return GenomicRangesList(granges, sequence_file_path)


def search_name(pattern, data):
    """Searches pattern in the data.

    Args:
        pattern (str): a regex pattern, r'string.
        data (str): a data to search pattern in.

    Returns:
        str or None: search result or None
            (if data is None or no entry was found).
    """
    if data is None:
        return None
    result = pattern.search(data)
    if result is None:
        return None
    return result[1]


check_name_regex = re.compile(".*\(.*\).*")


def parse_gtf_start(field):
    """Parses GTF start field.

    Args:
        field (str): gtf start field.

    Returns:
        int: start value.
    """
    return int(field) - 1


def parse_gtf_end(field):
    """Parses GTF end field.

    Args:
        field (str): gtf end field.

    Returns:
        int: end value.
    """
    return int(field)


def parse_phase(field):
    """Parses GTF/GFF phase field.

    Args:
        field (str): phase field.

    Returns:
        str or int: '.' or 0, 1, 2.

    Raises:
        IncorrectPhase: in case phase is not
            one of '.', 0, 1, 2.
    """
    if field == '.':
        return field
    elif field in ('0', '1', '2'):
        return int(field)
    raise IncorrectPhase(field)


def parse_gtf_attributes(field):
    """Parses GTF attributes field.

    Args:
        field (str): attributes field.

    Returns:
        dict: parse attributes field in a
            key-value format.

    Raises:
        IncorrectGTFAttrs: in case of empty field or
            the field cannot be parsed due to incorrect format.
    """
    try:
        data = {tag: value.strip('"')
                for tag, value in (item.strip().split(' ', 1)
                                   for item in field.split(';')[:-1])}
        if len(data) == 0:
            raise IncorrectGTFAttrs(field)
        return data
    except Exception:
        raise IncorrectGTFAttrs(field)


gtf_fields = {'chrom': str,
              'source': str,
              'method': str,
              'start': parse_gtf_start,
              'end': parse_gtf_end,
              'score': parse_score,
              'strand': parse_strand,
              'phase': parse_phase,
              'attributes': parse_gtf_attributes}


gtf_fields_fast = {'chrom': str,
                   'source': str,
                   'method': str,
                   'start': parse_gtf_start,
                   'end': parse_gtf_end,
                   'score': parse_score,
                   'strand': parse_strand,
                   'phase': parse_phase,
                   'attributes': str}

gtf_comment = '#'


def gtf_parser(fileobj, verbose=False, sequence_file_path=None,
               name_regex=None, name_tag=None, parse_attributes=False):
    granges = list()
    if parse_attributes:
        colnames = gtf_fields.keys()
        dtypes = gtf_fields.values()
    else:
        colnames = gtf_fields_fast.keys()
        dtypes = gtf_fields_fast.values()
    comment = gtf_comment

    if name_regex is None:
        get_name = lambda i: i
    elif check_name_regex.match(name_regex) is None:
        raise ValueError
    else:
        pattern = re.compile(name_regex)
        get_name = lambda i: search_name(pattern, i)

    if parse_attributes:
        extract_name = lambda i: get_name(i.get(name_tag))
    else:
        extract_name = get_name

    for line_no, line in tqdm(enumerate(fileobj), disable=not verbose):
        if line.startswith(comment):
            continue
        try:
            data = tuple(line.strip().split('\t'))
            if len(data) != 9:
                raise InconsistentGTF(line)
            record = {key: func(item)
                      for key, func, item in zip(colnames, dtypes, data)}
            record['name'] = extract_name(record['attributes'])
            granges.append(GenomicRange(**record))

        except Exception:
            raise IncorrectLine(line, line_no + 1)

    return GenomicRangesList(granges, sequence_file_path)


def parse_gff_start(field):
    return int(field) - 1


def parse_gff_end(field):
    return int(field)


def parse_gff_attributes(field):
    try:
        data = {tag: value
                for tag, value in (item.split('=')
                                   for item in field.split(';'))}
        if len(data) == 0:
            raise IncorrectGFFAttrs(field)
        return data
    except Exception:
        raise IncorrectGFFAttrs(field)


gff_fields = {'chrom': str,
              'source': str,
              'type': str,
              'start': parse_gff_start,
              'end': parse_gff_end,
              'score': parse_score,
              'strand': parse_strand,
              'phase': parse_phase,
              'attributes': parse_gff_attributes}

gff_fields_fast = {'chrom': str,
                   'source': str,
                   'type': str,
                   'start': parse_gff_start,
                   'end': parse_gff_end,
                   'score': parse_score,
                   'strand': parse_strand,
                   'phase': parse_phase,
                   'attributes': str}

gff_comment = "#"


def gff_parser(fileobj, verbose=False, sequence_file_path=None,
               name_regex=None, name_tag=None, parse_attributes=False):
    granges = list()
    if parse_attributes:
        colnames = gff_fields.keys()
        dtypes = gff_fields.values()
    else:
        colnames = gff_fields_fast.keys()
        dtypes = gff_fields_fast.values()
    comment = gff_comment

    if name_regex is None:
        get_name = lambda i: i
    elif check_name_regex.match(name_regex) is None:
        raise ValueError
    else:
        pattern = re.compile(name_regex)
        get_name = lambda i: search_name(pattern, i)

    if parse_attributes:
        extract_name = lambda i: get_name(i.get(name_tag))
    else:
        extract_name = get_name

    for line_no, line in tqdm(enumerate(fileobj), disable=not verbose):
        if line.startswith(comment):
            continue
        try:
            data = tuple(line.strip().split('\t'))
            if len(data) != 9:
                raise InconsistentGFF(line)
            record = {key: func(item)
                      for key, func, item in zip(colnames, dtypes, data)}
            record['name'] = extract_name(record['attributes'])
            granges.append(GenomicRange(**record))

        except Exception:
            raise IncorrectLine(line, line_no + 1)
    return GenomicRangesList(granges, sequence_file_path)


def check_bed3(line):
    try:
        record = line.strip().split('\t')
        if len(record) != 3:
            raise InconsistentBED3(line)
        dummy = {key: func(item)
                 for (key, func), item in zip(bed_fields.items(), record)}
        return True
    except Exception:
        return False


def check_bed6(line):
    try:
        record = line.strip().split('\t')
        if len(record) != 6:
            raise InconsistentBED6(line)
        dummy = {key: func(item)
                 for (key, func), item in zip(bed_fields.items(), record)}
        return True
    except Exception:
        return False


def check_bed12(line):
    try:
        record = line.strip().split('\t')
        if len(record) != 12:
            raise InconsistentBED12(line)
        dummy = {key: func(item)
                 for (key, func), item in zip(bed_fields.items(), record)}
        return True
    except Exception:
        return False


def check_gtf(line):
    try:
        record = line.strip().split('\t')
        if len(record) != 9:
            raise InconsistentGTF(line)
        dummy = {key: func(item)
                 for (key, func), item in zip(gtf_fields.items(), record)}
        return True
    except Exception:
        return False


def check_gff(line):
    try:
        record = line.strip().split('\t')
        if len(record) != 9:
            raise InconsistentGFF(line)
        dummy = {key: func(item)
                 for (key, func), item in zip(gff_fields.items(), record)}
        return True
    except Exception:
        return False


def annotation_sniffer(fileobj, comment='#', separator='\t'):
    line = fileobj.readline()

    while line.startswith(comment):
        line = fileobj.readline()
    if line == '':
        raise EmptyAnnotation
    record = line.strip().split(separator)

    if len(record) == 3:
        if check_bed3(line):
            fileobj.seek(0)
            return 'bed3'
        else:
            raise UnrecognizedFormat(3)
    elif len(record) == 6:
        if check_bed6(line):
            fileobj.seek(0)
            return 'bed6'
        else:
            raise UnrecognizedFormat(6)
    elif len(record) == 9:
        if check_gtf(line):
            fileobj.seek(0)
            return 'gtf'
        elif check_gff(line):
            fileobj.seek(0)
            return 'gff'
        else:
            raise UnrecognizedFormat(9)
    elif len(record) == 12:
        if check_bed12(line):
            fileobj.seek(0)
            return 'bed12'
        else:
            raise UnrecognizedFormat(12)
    else:
        raise UnrecognizedFormat(len(record))


def parse_annotation(fileobj, verbose=False, sequence_file_path=None,
                     name_regex=None, name_tag=None, parse_attributes=False):
    annotation_type = annotation_sniffer(fileobj)
    if annotation_type == 'bed3':
        return bed3_parser(fileobj, verbose, sequence_file_path)
    elif annotation_type == 'bed6':
        return bed6_parser(fileobj, verbose, sequence_file_path)
    elif annotation_type == 'bed12':
        return bed12_parser(fileobj, verbose, sequence_file_path)
    elif annotation_type == 'gtf':
        return gtf_parser(fileobj, verbose, sequence_file_path,
                          name_regex, name_tag, parse_attributes)
    elif annotation_type == 'gff':
        return gff_parser(fileobj, verbose, sequence_file_path,
                          name_regex, name_tag, parse_attributes)
