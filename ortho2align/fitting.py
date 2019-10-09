import abc
import numpy as np
from scipy.stats import rv_histogram, gaussian_kde


class AbstractFitter(metaclass=abc.ABCMeta):
    """Abstract class for empirical distribution estimation."""

    @property
    @abc.abstractmethod
    def estimator(self):
        """Abstract property for estimator object."""
        pass

    @abc.abstractmethod
    def pdf(self, data):
        """Returns probability density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                pdf with fitted empirical distribution.
        Returns:
            (number or np.array) pdf values of the items in data.
        """
        pass

    @abc.abstractmethod
    def cdf(self, data):
        """Returns cumulative density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                cdf with fitted empirical distribution.
        Returns:
            (number or np.array) cdf values of the items in data.
        """
        pass

    @abc.abstractmethod
    def sf(self, data):
        """Returns survival function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                sf with fitted empirical distribution.
        Returns:
            (number or np.array) sf values of the items in data.
        """
        pass

    @abc.abstractmethod
    def ppf(self, data):
        """Returns percent point function of the items in data.
        
        Args:
            data (number or iterable): data to estiamate
                ppf with fitted empirical distribution.
        Returns:
            (number or np.array) ppf values of the items in data.
        """
        pass

    @abc.abstractmethod
    def isf(self, data):
        """Returns inverse survival function of the items in data.
        
        Args:
            data (number or iterable): data to estiamate
                isf with fitted empirical distribution.
        Returns:
            (number or np.array) isf values of the items in data.
        """
        pass

    @classmethod
    def __subclasshook__(cls, C):
        """Magic method for issubclass() function."""
        if cls is AbstractFitter:
            attrs = set(dir(C))
            if set(cls.__abstractmethods__) <= attrs:
                return True

        return NotImplemented


def inverse_approximator(value, function, left, right, epsilon, increasing=True, iterations=100):
    """Finds argument of monotonic function given its result.
    
    Utilizes binary search algorithm to find the best argument
    so that the function result of it is the closest to the
    providied value. Search range is provided; function must
    be monotonic on this search range.
    
    Args:
        value (int, float): result of function.
        function (callable): function to calculate value.
        left (int, float): left boundary of the search range.
        right (int, float): right boundary of the search range.
        epsilon (int, float): the maximum difference between value
            and calculated result of function during the search to
            stop it.
        increasing (bool): whether the function is increasing on the
            search range or not (default: True).
        iterations (int): number of search iterations to perform
            until halt (default: 100).
    
    Returns:
        (int, float) the best argument found.
    
    """
    if right <= left:
        raise ValueError(f"`right` must be greater than `left`.")
    point = (right + left) / 2
    result = function(point)
    counter = 0
    while abs(result - value) > epsilon and counter < iterations:
        if result == value:
            return point
        elif increasing == (result > value):
            right = point
            point = (right + left) / 2
            result = function(point)
        else:
            left = point
            point = (right + left) / 2
            result = function(point)
        counter += 1
    return point


class HistogramFitter(AbstractFitter):
    """Fits empirical data with histogram and produces continuous distribution.
    
    Attributes:
        data (collection): numerical data used to fit distribution.
        estimator (scipy.rv_histogram) histogram instance used to calculate
            statistics.
    """

    def __init__(self, data):
        """Initializes a HistogramFitter instance.

        Args:
            data (iterable): empirical data to be fitted.

        Returns:
            None.
        """
        self.data = data
        self.hist = np.histogram(self.data,
                                 bins=len(self.data) // 2,
                                 density=True)
        self._estimator = rv_histogram(self.hist)

    @property
    def estimator(self):
        """Returns _estimator property."""
        return self._estimator

    def pdf(self, data):
        """Returns probability density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                pdf with fitted empirical distribution.
        Returns:
            (number or np.array) pdf values of the items in data.
        """
        return self.estimator.pdf(data)

    def cdf(self, data):
        """Returns cumulative density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                cdf with fitted empirical distribution.
        Returns:
            (number or np.array) cdf values of the items in data.
        """
        return self.estimator.cdf(data)

    def sf(self, data):
        """Returns survival function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                sf with fitted empirical distribution.
        Returns:
            (number or np.array) sf values of the items in data.
        """
        return self.estimator.sf(data)

    def ppf(self, data):
        """Returns percent point function of the items in data.
        
        Args:
            data (number or iterable): data to estiamate
                ppf with fitted empirical distribution.
        Returns:
            (number or np.array) ppf values of the items in data.
        """
        return self.estimator.ppf(data)

    def isf(self, data):
        """Returns inverse survival function of the items in data.
        
        Args:
            data (number or iterable): data to estiamate
                isf with fitted empirical distribution.
        Returns:
            (number or np.array) isf values of the items in data.
        """
        return self.estimator.isf(data)


class KernelFitter(AbstractFitter):
    """Fits empirical data with KDE and produces continuous distribution.
    
    Attributes:
        data (collection): numerical data used to fit distribution.
        estimator (scipy.gaussian_kde) KDE instance used to calculate
            statistics.
    """

    def __init__(self, data):
        """Initializes KernelFitter instance.
        
        Args:
            data (collection): numerical data to be fitted.
        
        Returns:
            None
        """
        self.data = data
        self._estimator = gaussian_kde(self.data)

    @property
    def estimator(self):
        """Returns _estimator property."""
        return self._estimator

    def pdf(self, data):
        """Returns probability density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                pdf with fitted empirical distribution.
        Returns:
            (number or np.array) pdf values of the items in data.
        """
        return self.estimator.pdf(data)

    def cdf(self, data):
        """Returns cumulative density function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                cdf with fitted empirical distribution.
        Returns:
            (number or np.array) cdf values of the items in data.
        """
        if isinstance(data, int) or isinstance(data, float):
            return self.estimator.integrate_box_1d(np.NINF, data)
        return np.array([self.estimator.integrate_box_1d(np.NINF, x)
                         for x in data])

    def sf(self, data):
        """Returns survival function of the items in data.

        Args:
            data (number or iterable): data to estiamate
                sf with fitted empirical distribution.
        Returns:
            (number or np.array) sf values of the items in data.
        """
        return 1 - self.cdf(data)

    def ppf(self, data, epsilon='auto'):
        """Returns percent point function of the items in data.
        
        Utilizes binary search to find ppf of values in data.
        
        Args:
            data (number or iterable): data to estiamate
                ppf with fitted empirical distribution.
            epsilon (int, float, str): the maximum difference between data and ppf(points)
                during the search to stop it. If 'auto', than selects average size
                between data points used to fit distribution (default: 'auto').
        Returns:
            (number or np.array) ppf values of the items in data.
        """
        if isinstance(data, int) or isinstance(data, float):
            data = [data]
        points = list()
        if epsilon == 'auto':
            epsilon = (max(self.data) - min(self.data)) / self.estimator.n
        elif not (isinstance(epsilon, int) or isinstance(epsilon, float)):
            raise ValueError('Provided epsilon is not an "auto", int or float instance.')
        for p in data:
            point = inverse_approximator(p,
                                         self.cdf,
                                         min(self.data),
                                         max(self.data),
                                         epsilon)
            points.append(point)
        if len(points) == 1:
            return points[0]
        return np.array(points)

    def isf(self, data, epsilon='auto'):
        """Returns inverse survival function of the items in data.
        
        Utilizes binary search to find isf of values in data.
        
        Args:
            data (number or iterable): data to estiamate
                isf with fitted empirical distribution.
            epsilon (int, float, str): the maximum difference between data and isf(points)
                during the search to stop it. If 'auto', than selects average size
                between data points used to fit distribution (default: 'auto').
        Returns:
            (number or np.array) isf values of the items in data.
        """
        if isinstance(data, int) or isinstance(data, float):
            data = [data]
        points = list()
        if epsilon == 'auto':
            epsilon = (max(self.data) - min(self.data)) / self.estimator.n
        elif not (isinstance(epsilon, int) or isinstance(epsilon, float)):
            raise ValueError('Provided epsilon is not an "auto", int or float instance.')
        for p in data:
            point = inverse_approximator(p,
                                         self.sf,
                                         min(self.data),
                                         max(self.data),
                                         epsilon,
                                         increasing=False)
            points.append(point)
        if len(points) == 1:
            return points[0]
        return np.array(points)
