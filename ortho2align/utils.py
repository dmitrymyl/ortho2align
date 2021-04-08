from itertools import groupby
from collections import Counter


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


def str_int_to_split(str_int):
    """Produces split representation of str int for slplot.

    Splits every int in two parts according to the following
    scheme:
    * one-digit numbers: (0, number);
    * two-digit numbers: (the first digit, the second digit);
    * three-digit numbers: (the first digit, remaining digits);
    * four-digit numbers and longer: (first digits, the last two digits).

    Args:
        str_int (str of int): str representation of an int.

    Returns:
        (tuple) 3 values: length of str_int, first digits and
            remaining digits.
    """
    if len(str_int) == 1:
        return 2, 0, str_int
    elif 2 <= len(str_int) <= 3:
        return len(str_int), str_int[0], str_int[1:]
    else:
        return len(str_int), str_int[:-2], str_int[-2:]


def slplot(data):
    """Returns a stem-and-leaf plot for given data.

    Args:
        data (sequence of ints): integer data to
            plot stem-and-leaf plot for.

    Returns:
        (str) stem-and-leaf plot of given data.
    """
    str_split_data = (str_int_to_split(str_int)
                      for str_int in (str(int(item))
                                      for item in sorted(data)))
    result = ((len_key, f"{digit_key}| {','.join(tuple(item[2] for item in digit_group))},")
              for len_key, len_group in groupby(str_split_data, lambda item: item[0])
              for digit_key, digit_group in groupby(len_group, lambda item: item[1]))
    str_result = '\n-----\n'.join('\n'.join(item[1] for item in group) for key, group in groupby(result, lambda item: item[0]))
    return str_result


def simple_hist(data):
    """Returns a simple histogram for the data."""
    counts = Counter(data)
    result = {key: value
              for key, value in sorted(counts.items(),
                                       key=lambda i: (i[1], i[0]))}
    str_result = "\n".join([f"{key}: {value}" for key, value in result.items()])
    return str_result
