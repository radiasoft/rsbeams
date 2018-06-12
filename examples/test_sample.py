import matplotlib.pyplot as plt
from rsbeams.rsstats.sample import sample, calculate_cdf
import numpy as np

"""
Testing rsbeams.rsstats.sample using a Gaussian distribution.

The sample function of the sample module can be used to draw a random set of points
from an arbitrary distribution given by some function.
"""


# Fixed seed for testing
np.random.seed(12345)


# Function to draw sample from (Gaussian Distribution)
def gaussian(x, sigma):
    return np.exp(-x**2 / (2. * sigma**2)) / np.sqrt(2. * np.pi * sigma**2)

# Values for sampling
a, b = -1., 1.
gaussian_width = 0.09
test_values = sample(10000, gaussian, a, b, args=(gaussian_width,))

# Plot result
gauss, ax1 = plt.subplots(1, 1, figsize=(10, 10))

# Histogram of sampled values
vals, bins = np.histogram(test_values, bins=128, range=(a, b))
vals = np.array(vals) / np.max(vals).astype(float)
centers = []
for i in range(bins.size - 1):
    centers.append((bins[i + 1] + bins[i]) / 2.)
ax1.plot(centers, vals)

# Plot of actual distribution function
xval = np.linspace(a, b, 10000)
gauss_max = np.max(gaussian(xval, gaussian_width))
ax1.plot(xval, gaussian(xval, gaussian_width) / gauss_max)

plt.show()
