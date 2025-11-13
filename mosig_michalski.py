# -*- coding: utf-8 -*-
"""
Uso do algoritmo de Mosig-Michalski para computar integrais de Sommerfeld.

Dada a seguinte integral que considera uma corrente em uma região circular de
raio $\tau$, onde $z$ e $z'$ são distâncias perpendiculares à estratificação,
do campo e da fonte, respectivamente, e $\rho$ é a distância lateral entre esses
pontos.

# Math: \int_0^\infty \tilde{G}(k_\rho; z|z') \, J_\nu(k_\rho \rho) \, J_\sigma(k_\rho \tau) \, k_\rho \, \mathrm{d}k_\rho

Estas integrais são frequentemente oscilatórias e divergentes com polos
(singularidades) no eixo real em $a_1$, $a_2$, etc. Para computar a cauda desta
integral mais rapidamente, faz-se partição e extrapolação (via algoritmo de
Mosig-Michalski) da série infinita que aproxima a integral.

# Referências
- MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
- MICHALSKI, K. A. On the efficient evaluation of integral arising in the sommerfeld halfspace problem. In: IEE Proceedings H (Microwaves, Antennas and Propagation). IET Digital Library, 1985. p. 312-318.
"""

import numpy as np
from scipy.integrate import tanhsinh


def part_extrap(fun, a, q, z, alpha, tol=1e-6, kmax=10, u=1):
    """
    Numerical integration of Sommerfeld integral tails with Partition-Extrapolation
    via the Mosig–Michalski method.

    The expected form of the integral is the following, where ``tau`` is the
    radius of the circular region that contains nonzero current, ``z`` and ``z'``
    are the perpendicular distances to the medium stratification of the field
    and the source, respectively, and ``rho`` is the lateral distance between
    those points:

    .. math:: \\int_0^\\infty \\tilde{G}(k_\\rho; z|z') \\, J_\\nu(k_\\rho \\rho) \\, J_\\sigma(k_\\rho \\tau) \\, k_\\rho \\, \\mathrm{d}k_\\rho

    Parameters
    ----------
    fun : callable
        The integrand function f(x) to be integrated over each partition subinterval.
    a : float
        The lower limit of integration (typically the start of the tail, after any head subtraction).
    q : float
        Partition step size; the width of each subinterval (typically matched to Bessel oscillation period or chosen for convergence).
    z : float
        Tail asymptotic decay parameter (often related to geometry or layer separation; used in remainder estimate).
    alpha : float
        Decay exponent in analytic remainder estimate for the extrapolation transformation.
    tol : float, optional
        Desired relative tolerance for extrapolation stopping criterion.
    kmax : int, optional
        Maximum number of partitions/subintervals for extrapolation.
    u : int, optional
        Transformation parameter controlling the extrapolation's sensitivity (1 for logarithmic convergence, 2 otherwise).

    Returns
    -------
    val : complex
        Accelerated estimate of the semi-infinite tail integral.
    error_estimate : float
        Estimate of the extrapolation absolute error.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.integrate import tanhsinh
    >>> from scipy.special import jv, gamma
    >>> v = 2
    >>> p = 1
    >>> a = 5.13562
    >>> q = np.pi / p
    >>> z = 0
    >>> tol = 1e-9
    >>> alpha = 1 / 2 - v
    >>> kmax = 10
    >>> u = 2
    >>> def fun(x):
    ...    return np.exp(-x * z) * jv(v, x * p) * x**v
    >>> y1 = ((2 * p) ** v* gamma(v + 1 / 2) / (np.sqrt(np.pi) * (z**2 + p**2) ** (v + 1 / 2)))
    >>> res = tanhsinh(fun, 0, a)
    >>> assert res.success
    >>> y2 = res.integral
    >>> Sn = y1 - y2  # "exact" value of the integral tail
    >>> val, E = part_extrap(fun, a, q, z, alpha, tol, kmax, u)
    >>> assert abs(val - Sn) < tol
    """
    X = np.zeros(kmax + 2)
    X[0] = a
    partial_sums = np.zeros(kmax + 2, dtype=complex)
    partial_errors = np.zeros_like(partial_sums)

    for k in range(kmax + 1):
        X[k + 1] = X[k] + q
        res = tanhsinh(fun, X[k], X[k + 1])
        if res.success is False:
            print(f"Numerical integration was not successful at partition {k}")
            print(f"tanhsinh status: {res.status}")

        partial_sums[k + 1] = res.integral + partial_sums[k]
        partial_errors[k + 1] = res.error + partial_errors[k]

    R = np.array(partial_sums, copy=True)
    old0 = 0.0
    old1 = 0.0
    val = 0.0
    error_estimate = 0.0
    for k in range(kmax + 1):
        if k > 0:
            # analytical remainder estimates
            Gk = -np.exp(-q * z) * (X[k] / X[k + 1]) ** alpha
        else:
            Gk = 0.0

        for j in range(1, k + 1):
            d = X[k - j + 2] - X[k - j + 1]
            eta = Gk / (1 + u * (j - 1) * d / X[k - j + 1])
            R[k - j + 1] = (R[k - j + 2] - eta * R[k - j + 1]) / (1 - eta)

        val = R[1]
        error_estimate = max(abs(val - old0), abs(val - old1))
        if k > 1 and error_estimate < tol * abs(val):
            break
        old0 = old1
        old1 = val

    return val, error_estimate


if __name__ == "__main__":
    from scipy.special import jv, gamma

    v = 2
    p = 1
    a = 5.13562
    q = np.pi / p
    z = 0
    tol = 1e-9
    alpha = 1 / 2 - v
    kmax = 10
    u = 2

    # Define the integrand function
    def fun(x):
        return np.exp(-x * z) * jv(v, x * p) * x**v

    # Reference ("exact") value
    y1 = (
        (2 * p) ** v
        * gamma(v + 1 / 2)
        / (np.sqrt(np.pi) * (z**2 + p**2) ** (v + 1 / 2))
    )
    res = tanhsinh(fun, 0, a)
    assert res.success
    y2 = res.integral
    Sn = y1 - y2  # "exact" value of the integral tail

    val, E = part_extrap(fun, a, q, z, alpha, tol, kmax, u)

    assert abs(val - Sn) < tol
