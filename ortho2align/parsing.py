import re
from tqdm.auto import tqdm
from .genomicranges import GenomicRange, GenomicRangesList


class ParserException(Exception):
    pass


class IncorrectStrand(ParserException):
    def __init__(self, strand):
        self.strand = strand

    def __str__(self):
        return f"Incorrect strand value: {self.strand}."


class IncorrectItemRgb(ParserException):
    def __init__(self, itemRgb):
        self.itemRgb = itemRgb

    def __str__(self):
        return f"Incorrenct itemRgb value: {self.itemRgb}"


class IncorrectBlockSizes(ParserException):
    def __init__(self, blockSizes):
        self.blockSizes = blockSizes

    def __str__(self):
        return f"Incorrect blockSizes value: {self.blockSizes}"


class InconsistentBlockSizes(IncorrectBlockSizes):
    def __init__(self, blockSizes, blockCount):
        self.blockSizes = blockSizes
        self.blockCount = blockCount

    def __str__(self):
        return f"blockSizes field '{self.blockSizes}' is inconsistent with blockCount field '{self.blockCount}'."


class IncorrectBlockStarts(ParserException):
    def __init__(self, blockStarts):
        self.blockStarts = blockStarts

    def __str__(self):
        return f"Incorrect blockStarts value: {self.blockStarts}"


class InconsistentBlockStarts(IncorrectBlockStarts):
    def __init__(self, blockStarts, blockCount):
        self.blockStarts = blockStarts
        self.blockCount = blockCount

    def __str__(self):
        return f"blockStarts field '{self.blockStarts}' is inconsistent with blockCount field '{self.blockCount}'."


class IncorrectPhase(ParserException):
    def __init__(self, phase):
        self.phase = phase

    def __str__(self):
        return f"Incorrect phase value: {self.phase}"


class IncorrectGTFAttrs(ParserException):
    def __init__(self, attrs):
        self.attrs = attrs

    def __str__(self):
        return f"Incorrect GTF attributes: {self.attrs}"


class IncorrectGFFAttrs(ParserException):
    def __init__(self, attrs):
        self.attrs = attrs

    def __str__(self):
        return f"Incorrect GFF attributes: {self.attrs}"


class InconsistentBED3(ParserException):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED3 file format in the following line: {self.line}"


class InconsistentBED6(ParserException):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED6 file format in the following line: {self.line}"


class InconsistentBED12(ParserException):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with BED12 file format in the following line: {self.line}"


class InconsistentGTF(ParserException):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with GTF file format in the following line: {self.line}"


class InconsistentGFF(ParserException):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return f"Number of fields is inconsistent with GFF file format in the following line: {self.line}"


class IncorrectLine(ParserException):
    def __init__(self, line, line_no):
        self.line = line
        self.line_no = line_no

    def __str__(self):
        return f"Line No. {self.line_no} is incorrect for the selected file format: {self.line}"


class EmptyAnnotation(ParserException):
    def __str__(self):
        return "Provided annotation file is empty or consist only of comment lines."


class UnrecognizedFormat(ParserException):
    def __init__(self, fields):
        self.fields = fields

    def __str__(self):
        return f"Unrecognized annotation format with {self.fields} fields."


def parse_bed_start(field):
    return int(field)


def parse_bed_end(field):
    return int(field)


def parse_score(field):
    try:
        a = float(field)
        if a.is_integer():
            return int(a)
        else:
            return a
    except ValueError:
        return field


def parse_strand(field):
    if field in ('.', '+', '-'):
        return field
    raise IncorrectStrand(field)


def parse_itemRgb(field):
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
    try:
        record = tuple(int(i) for i in field.split(',') if i != '')
        return record
    except Exception:
        raise IncorrectBlockSizes(field)


def parse_blockStarts(field):
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
                raise IncorrectBlockStarts(record['blockStarts'], record['blockCount'])
            granges.append(GenomicRange(**record))
        except Exception:
            raise IncorrectLine(line, line_no + 1)

    return GenomicRangesList(granges, sequence_file_path)


def search_name(pattern, data):
    if data is None:
        return None
    result = pattern.search(data)
    if result is None:
        return None
    return result[1]


check_name_regex = re.compile(".*\(.*\).*")


def parse_gtf_start(field):
    return int(field) - 1


def parse_gtf_end(field):
    return int(field)


def parse_phase(field):
    if field == '.':
        return field
    elif field in ('0', '1', '2'):
        return int(field)
    raise IncorrectPhase(field)


def parse_gtf_attributes(field):
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
