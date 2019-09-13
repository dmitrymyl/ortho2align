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
        pass

    @abc.abstractmethod
    def isf(self, data):
        pass

    @classmethod
    def __subclasshook__(cls, C):
        """Magic method for issubclass() function."""
        if cls is AbstractFitter:
            attrs = set(dir(C))
            if set(cls.__abstractmethods__) <= attrs:
                return True

        return NotImplemented


class HistogramFitter(AbstractFitter):
    """Fits empirical data with histogram and produces continuous distribution."""

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
        return self.estimator.ppf(data)

    def isf(self, data):
        return self.estimator.isf(data)


class KernelFitter(AbstractFitter):

    def __init__(self, data):
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
        # if isinstance(data, int) or isinstance(data, float):
        #     return self.estimator.integrate_box_1d(data, np.Inf)
        # return np.array([self.estimator.integrate_box_1d(x, np.Inf)
        #                  for x in data])
        return 1 - self.cdf(data)

    def ppf(self, data):
        if isinstance(data, int) or isinstance(data, float):
            data = [data]
        points = list()
        for p in data:
            left = min(data)
            right = max(data)
            point = (right + left) / 2
            epsilon = (left - right) / self.estimator.n
            result = self.cdf(point)
            while abs(result - p) > epsilon:
                if result > p:
                    right = point
                    point = (right + left) / 2
                    result = self.cdf(point)
                    epsilon = (left - right) / self.estimator.n  # TODO: refactor
                else:
                    left = point
                    point = (right + left) / 2
                    result = self.cdf(point)
                    epsilon = (left - right) / self.estimator.n  # TODO: refactor
            points.append(point)
        if len(points) == 1:
            return points[0]
        return np.array(points)

    def isf(self, data):
        return 1 - self.ppf(data)
