import numpy as np
from pathos.multiprocessing import Pool, cpu_count
from scipy.integrate import quad
from scipy.special import fresnel


class CentroidPosition:
    def __init__(self, N, Z, nu0, mu):
        self.N = N
        self.Z = Z
        self.nu0 = nu0
        self.mu = mu

    def _reduced_integrand(self, a, n):
        order = 1
        advance = 0
        for m in self.mu:
            advance += m * a ** order / (2. * np.pi * n) ** (order - 1)
            order += 1

        coeff = self.Z / (2 * n)
        const_slip = 2 * np.pi * self.nu0 * n
        angular_term = np.cos(const_slip) * np.cos(advance) + np.sin(const_slip) * np.sin(advance)

        # Calculate cutoff if a is float or array
        try:
            maxa = 1. * 2 * np.pi * n
            if a <= maxa:
                distr = angular_term / 1. / np.pi
            else:
                distr = 0.
        except ValueError:
            maxa = np.ones_like(a, dtype='float') * 2 * np.pi * n
            distr = angular_term / 1. / np.pi * np.less(a, maxa)

        return coeff * distr

    def integrate_any_order(self, turn=None):
        if turn is not None:
            n = turn
        else:
            n = self.N
        if n == 0:
            return self.Z

        # noinspection PyTupleAssignmentBalance
        result, _ = quad(self._reduced_integrand,
                         0, 2 * np.pi * n,
                         args=n)
        return result

    def integrate_first_order(self, turn=None):
        if turn is not None:
            n = turn
        else:
            n = self.N
        if n == 0:
            return self.Z

        xN = self.Z / (2. * np.pi * n * self.mu[0]) * \
             (np.cos(2 * np.pi * self.nu0 * n) * np.sin(2 * np.pi * n * self.mu[0]) +
              2. * np.sin(2 * np.pi * self.nu0 * n) * np.sin(
                  np.pi * n * self.mu[0]) ** 2)

        return xN

    def integrate_second_order(self, turn=None):
        if turn is not None:
            n = turn
        else:
            n = self.N
        if n == 0:
            return self.Z

        def integrand(u, N):
            fS, fC = fresnel((self.mu[0] * N * np.pi + self.mu[1] * u) / np.sqrt(self.mu[1] * N * np.pi**2))
            term1 = np.cos(np.pi * self.mu[0]**2 * N / (2. * self.mu[1]) + 2. * np.pi * self.nu0 * N)
            term2 = np.sin(np.pi * self.mu[0]**2 * N / (2. * self.mu[1]) + 2. * np.pi * self.nu0 * N)

            return fC * term1 + fS * term2

        xN = integrand(2 * np.pi * n, n) - integrand(0, n)

        return xN * self.Z / np.sqrt(4. * self.mu[1] * n)

    def calculate_centroids(self, p=None):
        if p:
            pool_size = p
        else:
            pool_size = cpu_count()

        pool = Pool(pool_size)

        #  attempt to speed things up by spreading out difficult integration values at the end of range
        #  appeared to not work
        #     x = []
        #     for i in range(cpu_count()):
        #         x += range(N)[i::4]

        if len(self.mu) == 1:
            integration_function = self.integrate_first_order
        elif len(self.mu) == 2:
            integration_function = self.integrate_second_order
        else:
            integration_function = self.integrate_any_order

        x = range(self.N)
        results = pool.map(integration_function, x)
        pool.close()

        return results
