from numpy import trapz


def calc_precision(TP, FP):
    """Calculates precision from TP and FP.

    Args:
        TP (int, float): number of true positive cases.
        FP (int, float): number of false positive cases.

    Returns:
        (float): precision value.
    """
    return TP / (TP + FP)


def calc_recall(TP, FN):
    """Calculates recall (sensitivity) from TP and FN.

    Args:
        TP (int, float): number of true positive cases.
        FN (int, float): number of false negative cases.

    Returns:
        (float): recall value.
    """
    return TP / (TP + FN)


def calc_specificity(TN, FP):
    """Calculates specificity (sensitivity) from TN and FP.

    Args:
        TN (int, float): number of true negative cases.
        FP (int, float): number of false positive cases.

    Returns:
        (float): specificity value.
    """
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
    return (TP + TN) / (TP + TN + FP + FN)


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


def trace_orthologs(granges_list, found_relname='found', real_relname='real'):
    for grange in granges_list:
        for found_grange in grange.relations[found_relname]:
            found_grange.relations['trace'] = found_grange.find_neighbours(grange.relations[real_relname])
        for real_grange in grange.relations[found_relname]:
            real_grange.relations['trace'] = real_grange.find_neighbours(grange.relations[real_relname])


def calc_ortholog_metrics(granges_list, found_relname='found', real_relname='real'):
    metric_names = ['TP', 'FP', 'TN', 'FN', 'JI', 'OC', 'QF', 'FF', 'RF']
    metrics = {name: list() for name in metric_names}
    for grange in granges_list:
        tp, tn, fp, fn, ji, oc, qf, ff, rf = [0 for _ in range(9)]
        if len(grange.relations[found_relname]) == 0 and \
           len(grange.relations[real_relname]) == 0:
            tn = 1
        elif len(grange.relations[found_relname]) > 0 and \
             len(grange.relations[real_relname]) > 0:
            tp = sum([len(found_grange.relations['trace'])
                      for found_grange in grange.relations[found_relname]])
            fp = len([found_grange
                      for found_grange in grange.relations[found_relname]
                      if len(found_grange.relations['trace']) == 0])
            fn = len([real_grange
                      for real_grange in grange.relations[real_relname]
                      if len(real_grange.relations['trace']) == 0])
            ji = [found_grange.calc_JI(real_grange)
                  for found_grange in grange.relations[found_relname]
                  for real_grange in found_grange.relations['trace']]
            if len(ji) == 0:
                ji = 0
            oc = [found_grange.calc_OC(real_grange)
                  for found_grange in grange.relations[found_relname]
                  for real_grange in found_grange.relations['trace']]
            if len(oc) == 0:
                oc = 0
            qf = [grange.calc_fraction(found_grange)
                  for found_grange in grange.relations[found_relname]]
            if len(qf) == 0:
                qf = 0
            ff = [real_grange.calc_fraction(found_grange)
                  for real_grange in grange.relations[real_relname]
                  for found_grange in real_grange.relations['trace']]
            if len(ff) == 0:
                ff = 0
            rf = [found_grange.calc_fraction(real_grange)
                  for found_grange in grange.relations[found_relname]
                  for real_grange in found_grange.relations['trace']]
            if len(rf) == 0:
                rf = 0
        elif len(grange.relations[found_relname]) > 0:
            fp = len(grange.relations[found_relname])
        else:  # len(grange.relations[real_relname]) > 0
            fn = len(grange.relations[real_relname])
        grange_metrics = [tp, fp, tn, fn, ji, oc, qf, ff, rf]
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


def calc_ensemble_metric(ortholog_metrics, metric_name):
    """Sums a metric from calc_ortholog_metrics."""
    if metric_name in ortholog_metrics:
        return sum(unnest_number_list(ortholog_metrics[metric_name]))
    else:
        raise ValueError(f'metric_name {metric_name} is not one of {ortholog_metrics.keys()}.')
