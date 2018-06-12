import numpy as np


def calculate_cdf(weights):
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
    return cdf_array


def sample(n, distribution, a, b, args=None, wn=1000):
    """
    Draws a set of samples from a distribution by calculating the cumulative distribution function
    and using the inverse transform sampling method.

    Args:
        n: (int) Number of points to sample from `distribution`.
        distribution: Distribution function to be sampled. Must be normalized such
            that the area of the distribution is 1. Must be a callable function of form f(x, *args),
            returning the probability distribution as a function of x with option parameters.
        a: Lower bound to be used in calculating weights for sampling.
        b: Upper bound to be used in calculating weights for sampling.
        args: (tuple) Optional set of arguments to be passed to `distribution`.
        wn: Number of points to use in calculating weights from `distribution`.

    Returns:
        (ndarray) of `n` points drawn from `distribution`.
    """
    assert callable(distribution), "`distribution` must be a callable function."
    x = np.linspace(a, b, wn)

    if args:
        try:
            args[0]
        except TypeError:
            args = (args,)
        weights = distribution(x, *args)
    else:
        weights = distribution(x)

    cdf_vals = calculate_cdf(weights)
    y = np.random.random(n)
    samples = np.searchsorted(cdf_vals, y)

    return x[samples]
