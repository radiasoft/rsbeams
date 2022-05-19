import numpy as np
from scipy import stats
from scipy import interpolate


def _calculate_cdf(weights: np.ndarray) -> np.ndarray:
    """
    Calculate the cumulative distribution function given a discrete set of weights for a probability function.

    Args:
        weights: (Iterable) Set of values representing weights for some distribution function (Ideally normalized
            to unity area).

    Returns:
        (ndarray) Array of values for the discrete CDF function represented by `weights`, with len(`weights`).
    """
    total = np.sum(weights)
    cdf_array = []
    cumulative = 0
    for weight in weights:
        cumulative += weight
        cdf_array.append(cumulative / total)
    return np.array(cdf_array)


def cdf_on_grid(distribution: callable, grid: np.ndarray, args: tuple = None) -> np.ndarray:
    """Cumulative distribution on a grid.
    Calculates the cumulative distribution function (CDF) on a grid of points.

    Args:
        distribution: Distribution function to be sampled. Must be normalized such
            that the area of the distribution is 1. Must be a callable function of form f(x, *args),
            returning the probability distribution as a function of x with option parameters.
        grid: (ndarray) Grid of points to calculate the CDF on.
        args: (tuple) Optional set of arguments to be passed to `distribution`.

    Returns:
        (ndarray) of `n` points drawn from `distribution`.

    """

    if args:
        try:
            args[0]
        except TypeError:
            args = (args,)
        weights = distribution(grid, *args)
    else:
        weights = distribution(grid)

    cdf_vals = _calculate_cdf(weights)
    return cdf_vals


def sample(n: int, distribution: callable, a: int, b: int,
           args: tuple = None, wn: int = 10000, method: str = 'random', invert: str = 'numeric') -> np.ndarray:
    """Sample from a 1D distribution.
    Draws a set of samples from a distribution by calculating the cumulative distribution function
    and using the inverse transform sampling method.

    Args:
        n: (int) Number of points to sample from `distribution`.
        distribution: Distribution function to be sampled. Must be normalized such
            that the area of the distribution is 1. Must be a callable function of form f(x, *args),
            returning the probability distribution as a function of x with option parameters.
        a: (float) Lower bound to be used for sampling.
        b: (float) Upper bound to be used for sampling.
        args: (tuple) Optional set of arguments to be passed to `distribution`.
        wn: (int) [10000] Number of points on [a, b) used in calculating weights from `distribution`
            and the samples that will be returned.
        method: (str) Determines how samples will be drawn. Should be either 'random' or 'halton'.
            If 'random' samples will be drawn from a uniform distribution. For 'halton' sampling uses
            pseudo-random Halton sequences to deterministically generate the samples.
        invert: (str) How inversion of the cumulative density function (CDF) will be performed.
            must be either 'numeric' or 'interpolate'. If 'numeric' then the inverse CDF is approximated
            on the grid used to calculate weights. The result will be sensitive to the setting of `wn` and,
            if used with method='halton' will not preserve the Halton sequence's low-discrepancy properties.
            Using 'interpolate' will calculate the inverse CDF from interpolation. This method is normally
            preferable unless the interpolation fails.

    Returns:
        (ndarray) of `n` points drawn from `distribution`.
    """
    assert callable(distribution), "`distribution` must be a callable function."
    x = np.linspace(a, b, wn)

    cdf_vals = cdf_on_grid(distribution, x, args=args)

    if method == 'random':
        y = np.random.random(n)
    elif method == 'halton':
        sampler = stats.qmc.Halton(d=1, scramble=True)
        y = sampler.random(n=n)
    else:
        raise NotImplementedError(f'{method} is not supported for method')

    if invert == 'numeric':
        samples = np.searchsorted(cdf_vals, y)
        return x[samples]
    elif invert == 'interpolate':
        icdf = interpolate.interp1d(cdf_vals, x, kind='cubic', fill_value='extrapolate')
        samples = icdf(y)
        return samples
    else:
        raise NotImplementedError(f'{invert} is not supported for invert')
