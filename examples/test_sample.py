import matplotlib.pyplot as plt
from scipy.special import gamma
from rsbeams.rsstats.sample import sample, calculate_cdf
import numpy as np

"""
Testing rsbeams.rsstats.sample using a Gaussian distribution and Gamma distribution.

The sample function of the sample module can be used to draw a random set of points
from an arbitrary distribution given by some function.
"""


# Fixed seed for testing
np.random.seed(12345)

# Plot setup
gauss, [ax1, ax2] = plt.subplots(1, 2, figsize=(16, 10))


# Function to draw sample from
def gaussian_distr(x, sigma):
    return np.exp(-x**2 / (2. * sigma**2)) / np.sqrt(2. * np.pi * sigma**2)


def gamma_distr(x, shape, scale):
    return x ** (shape - 1) * ((np.exp(-x / scale)) / (gamma(shape) * scale ** shape))

# Values for sampling
a, b = -1., 1.
gaussian_width = 0.09
test_values = sample(10000, gaussian_distr, a, b, args=(gaussian_width,))

ax1.set_title("Gaussian Distribution")
ax1.set_xlabel("x")
ax1.set_ylabel("Counts")
# Histogram of sampled Gaussian values
vals, bins = np.histogram(test_values, bins=128, range=(a, b))
vals = np.array(vals) / np.max(vals).astype(float)
centers = []
for i in range(bins.size - 1):
    centers.append((bins[i + 1] + bins[i]) / 2.)
ax1.plot(centers, vals)

# Plot of actual Gaussian distribution function
xval = np.linspace(a, b, 10000)
gauss_max = np.max(gaussian_distr(xval, gaussian_width))
ax1.plot(xval, gaussian_distr(xval, gaussian_width) / gauss_max)

# Parameters for Gamma distribution
shape, scale = 2., 1.
a, b = 0., 18.
gamma_test_values = sample(10000, gamma_distr, a, b, args=(shape, scale))


ax2.set_title("Gamma Distribution")
ax2.set_xlabel("x")
ax2.set_ylabel("Counts")
# Histogram of sampled Gamma values
vals, bins = np.histogram(gamma_test_values, bins=128, range=(a, b))
vals = np.array(vals) / np.max(vals).astype(float)
centers = []
for i in range(bins.size - 1):
    centers.append((bins[i + 1] + bins[i]) / 2.)
ax2.plot(centers, vals)

# Plot of actual Gamma distribution function
xval = np.linspace(a, b, 10000)
gamma_max = np.max(gamma_distr(xval, shape, scale))
ax2.plot(xval, gamma_distr(xval, shape, scale) / gamma_max)


plt.show()
