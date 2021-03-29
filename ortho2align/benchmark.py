from .utils import numberize
from .genomicranges import BaseGenomicRange
from numpy import trapz


def calc_precision(TP, FP):
    """Calculates precision from TP and FP.

    Args:
        TP (int, float): number of true positive cases.
        FP (int, float): number of false positive cases.

    Returns:
        (float): precision value.
    """
    if TP == 0:
        if FP == 0:
            return 1
        else:
            return 0
    return TP / (TP + FP)


def calc_recall(TP, FN):
    """Calculates recall (sensitivity) from TP and FN.

    Args:
        TP (int, float): number of true positive cases.
        FN (int, float): number of false negative cases.

    Returns:
        (float): recall value.
    """
    if TP == 0:
        if FN == 0:
            return 1
        else:
            return 0
    return TP / (TP + FN)


def calc_specificity(TN, FP):
    """Calculates specificity (sensitivity) from TN and FP.

    Args:
        TN (int, float): number of true negative cases.
        FP (int, float): number of false positive cases.

    Returns:
        (float): specificity value.
    """
    if TN == 0:
        if FP == 0:
            return 1
        else:
            return 0
    return TN / (TN + FP)


def calc_accuracy(TP, TN, FP, FN):
    """Calculates accuracy from TP, TN, FP and FN.

    Args:
        TP (int, float): number of true positive cases.
        TN (int, float): number of true negative cases.
        FP (int, float): number of false positive cases.
        FN (int, float): number of false negative cases.

    Returns:
        (float): accuracy value.
    """
    try:
        return (TP + TN) / (TP + TN + FP + FN)
    except ZeroDivisionError:
        return 1


def calc_auc(x, y):
    """Calculates area under the curve using trapezoid rule.

    Args:
        x (iterable): x coordinates. Must be monotonic increasing.
        y (iterable): y coordinates.

    Returns:
        (float): AUC value.

    Raises:
        ValueError: in case lengths of x and y aren't equal.
    """
    if len(x) != len(y):
        raise ValueError(f"Lengths of x and y aren't equal: {len(x)} and {len(y)}.")
    return trapz(y, x)


def trace_orthologs(granges_list, found_subject_relname='found_subject', real_relname='real'):
    for grange in granges_list:
        for found_subject_grange in grange.relations[found_subject_relname]:
            found_subject_grange.relations['trace'] = found_subject_grange.find_neighbours(grange.relations[real_relname])
        for real_grange in grange.relations[real_relname]:
            real_grange.relations['trace'] = real_grange.find_neighbours(grange.relations[found_subject_relname])


def calc_ortholog_metrics(granges_list,
                          found_query_relname='found_query',
                          found_subject_relname='found_subject',
                          found_query_subject_relanme='ortho_link',
                          real_relname='real',
                          tp_mode='all'):
    metric_names = ['TP', 'FP', 'TN', 'FN', 'JI', 'OC', 'QF', 'FF', 'RF', 'CF']
    metrics = {name: list() for name in metric_names}
    for grange in granges_list:
        if len(grange.relations[found_subject_relname]) == 0 and \
           len(grange.relations[real_relname]) == 0:
            tn = [1]
            tp, fp, fn = [[0] for _ in range(3)]
            ji, oc, qf, ff, rf, cf = [[] for _ in range(6)]
        elif (len(grange.relations[found_subject_relname]) > 0
              and len(grange.relations[real_relname]) > 0):
            if tp_mode == 'all':
                tp = [sum([len(found_grange.relations['trace'])
                           for found_grange in grange.relations[found_subject_relname]])]
            elif tp_mode == 'single':
                tp_count = sum([len(found_grange.relations['trace'])
                                for found_grange in grange.relations[found_subject_relname]])
                tp = [1 if tp_count > 0 else 0]
            else:
                raise ValueError('tp_mode must be one of "all", "single".')
            fp = [len([found_grange
                       for found_grange in grange.relations[found_subject_relname]
                       if len(found_grange.relations['trace']) == 0])]
            fn = [len([real_grange
                       for real_grange in grange.relations[real_relname]
                       if len(real_grange.relations['trace']) == 0])]
            tn = [0]
            ji = [found_grange.calc_JI(real_grange)
                  for found_grange in grange.relations[found_subject_relname]
                  for real_grange in found_grange.relations['trace']]
            oc = [found_grange.calc_OC(real_grange)
                  for found_grange in grange.relations[found_subject_relname]
                  for real_grange in found_grange.relations['trace']]
            qf = [grange.calc_fraction(found_grange)
                  for found_grange in grange.relations[found_query_relname]]
            rf = [real_grange.calc_fraction(found_grange)
                  for real_grange in grange.relations[real_relname]
                  for found_grange in real_grange.relations['trace']]
            ff = [found_grange.calc_fraction(real_grange)
                  for found_grange in grange.relations[found_subject_relname]
                  for real_grange in found_grange.relations['trace']]
            cf = [grange.calc_fraction(found_query_grange) * found_subject_grange.calc_fraction(real_grange)
                  for found_query_grange in grange.relations[found_query_relname]
                  for found_subject_grange in found_query_grange.relations[found_query_subject_relanme]
                  for real_grange in found_subject_grange.relations['trace']]
        elif len(grange.relations[found_subject_relname]) > 0:
            fp = [len(grange.relations[found_subject_relname])]
            qf = [grange.calc_fraction(found_grange)
                  for found_grange in grange.relations[found_query_relname]]
            tp, tn, fn = [[0] for _ in range(3)]
            ji, oc, ff, rf, cf = [[] for _ in range(5)]
        else:  # len(grange.relations[real_relname]) > 0
            fn = [len(grange.relations[real_relname])]
            tp, tn, fp = [[0] for _ in range(3)]
            ji, oc, rf, ff = [[0] for _ in range(4)]
            qf, cf = [[] for _ in range(2)]
        grange_metrics = [tp, fp, tn, fn, ji, oc, qf, ff, rf, cf]
        for name, value in zip(metric_names, grange_metrics):
            metrics[name].append(value)
    return metrics


def unnest_number_list(collection):
    """Unnests a collection composed of numbers and lists."""
    unnested = list()
    for value in collection:
        if isinstance(value, list):
            unnested += value
        else:
            unnested.append(value)
    return unnested


def unnest_list_with_empties(collection):
    unnested = list()
    for value in collection:
        if isinstance(value, list):
            if value:
                unnested += value
            else:
                unnested.append(0)
        else:
            unnested.append(value)
    return unnested


def calc_ensemble_metric(ortholog_metrics, metric_name):
    """Sums a metric from calc_ortholog_metrics."""
    if metric_name in ortholog_metrics:
        return sum(unnest_number_list(ortholog_metrics[metric_name]))
    else:
        raise ValueError(f'metric_name {metric_name} is not one of {ortholog_metrics.keys()}.')


blastn_replace_dict = {'q. start': 'qstart',
                       'q. end': 'qend',
                       's. start': 'sstart',
                       's. end': 'send',
                       'query length': 'qlen',
                       'subject length': 'slen',
                       'query acc.ver': 'qchrom',
                       'subject acc.ver': 'schrom'}


def parse_blastn_results(file_object, query_grange):
    """Parse blastn search results for benchmarking purposes."""
    line = file_object.readline()
    if not line.startswith('#'):
        raise ValueError("Provided file_object is not in accepted file format.")

    for _ in range(100):
        line = file_object.readline()
        if line.startswith("# Fields:"):
            break
    else:
        return [], []

    fields = line.lstrip("# Fields: ").rstrip("\n").split(", ")

    fields = [blastn_replace_dict.get(item, item)
              for item in fields]
    hit_number = int(file_object.readline().split(" ")[1])
    query_HSPs = list()
    subject_HSPs = list()
    for i in range(hit_number):
        hsp = file_object.readline().strip().split("\t")
        hsp = [numberize(item) for item in hsp]
        hsp = dict(zip(fields, hsp))
        query_hsp = BaseGenomicRange(query_grange.chrom,
                                     hsp['qstart'],
                                     hsp['qend'],
                                     name=query_grange.name)
        subject_hsp = BaseGenomicRange(hsp['schrom'],
                                       hsp['sstart'],
                                       hsp['send'],
                                       name=query_grange.name)
        if query_grange.strand != '-':
            query_hsp.start = query_hsp.start + query_grange.start
            query_hsp.end = query_hsp.end + query_grange.start
            query_hsp.strand = '+'
        else:
            query_hsp.start = query_grange.end - query_hsp.start
            query_hsp.end = query_grange.end - query_hsp.end
            query_hsp.strand = '-'
        # subject sequence
        if subject_hsp.start <= subject_hsp.end:
            subject_hsp.strand = query_hsp.strand
        else:
            subject_hsp.start = subject_hsp.end
            subject_hsp.end = subject_hsp.start
            subject_hsp.strand = "+" if query_hsp.strand == "-" else "-"

        query_hsp.start -= 1
        subject_hsp.start -= 1
        query_HSPs.append(query_hsp)
        subject_HSPs.append(subject_hsp)
    return query_HSPs, subject_HSPs
