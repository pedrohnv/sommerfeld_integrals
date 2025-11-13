#=
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
=#

using DoubleExponentialFormulas
using SpecialFunctions


"""
Sequence acceleation of an infinite series by extrapolation via the Mosig–Michalski algorithm (t-transformation weighted averages).

# Arguments
- partial_sums: Vector of partial sums of the sequence to be accelerated, where the n-th entry should be the partial sum up to index n.
- u: Transformation parameter controlling the extrapolation's sensitivity (1 for logarithmic convergence, 2 otherwise).

# Returns
- extrapolated_sum: The extrapolated sum, providing an improved estimate for the series true value beyond available terms.

# Examples
```julia
    partial_sum(n::Int) = cumsum([series_term(i) for i = 0:(n - 1)])
    
    series_term(i::Int) = (-1)^i / sqrt(i + 1)
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, 0.604898643421630)

    series_term(i::Int) = (4/5)^(i+1) / (i + 1)
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, log(5))

    series_term(i::Int) = 1 / (i + 1)^2
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, pi^2 / 6; atol = 1e-2)
```

# References
- MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
"""
function mosig_michalski(partial_sums, u = 1)
    kmax = length(partial_sums) - 1

    beta = 1
    X = (0:kmax) .+ beta  # equidistant nodes

    # Transformed sequence values at each extrapolation order k
    R = deepcopy(partial_sums)

    # remainder estimates
    wk0 = 1  # w_{k-1}
    wk1 = 1  # w_{k}

    local extrapolated_sum
    for k = 0:(kmax)
        if k >= 1
            wk1 = partial_sums[k + 1] - partial_sums[k]
            Gk = wk1 / wk0
        else
            Gk = 0
        end
        for j = 1:k
            d = X[k - j + 2] - X[k - j + 1]
            eta = Gk / (1 + u * (j - 1) * d / X[k - j + 1])
            R[k - j + 1] = (R[k - j + 2] - eta * R[k - j + 1]) / (1 - eta)
        end
        extrapolated_sum = R[1]
        wk0 = wk1
    end
    return extrapolated_sum
end


"""
Numerical integragion of Sommerfeld integral tails with Partition-Extrapolation
via the Mosig-Michalski method.

The expected form of the integral is the following, where ``\\tau`` is the
radius of the circular region that contains nonzero current, ``z`` and ``z'``
are the perpendicular distances to the medium stratification of the field
and the source, respectively, and ``\\rho`` is the lateral distance between
those points.

```math
\\int_0^\\infty \\tilde{G}(k_\\rho; z|z') \\, J_\\nu(k_\\rho \\rho) \\, J_\\sigma(k_\\rho \\tau) \\, k_\\rho \\, \\mathrm{d}k_\\rho
```

# Parameters
- fun: The integrand function ``f(x)`` to be integrated over each partition subinterval.
- a: The lower limit of integration (typically the start of the tail, after any head subtraction).
- q: Partition step size; the width of each subinterval (typically matched to Bessel oscillation period or chosen for balance between cost and convergence).
- z: Tail asymptotic decay parameter (often related to physical geometry or stratification; used in remainder estimate).
- alpha: Decay exponent in analytic remainder estimate for the extrapolation transformation.
- tol: Desired relative tolerance for extrapolation stopping criterion.
- kmax: Maximum number of partitions/subintervals for extrapolation.
- u: Transformation parameter controlling the extrapolation's sensitivity (1 for logarithmic convergence, 2 otherwise).

# Returns
- val: Accelerated estimate of the semi-infinite tail integral.
- error_estimate: Estimate of extrapolation absolute error.

# Examples
```julia
    v = 2
    p = 1
    a = 5.13562
    q = pi / p
    z = 0
    tol = 1e-9
    alpha = 1/2 - v
    kmax = 10
    u = 2

    fun = x -> exp(-x * z) * besselj(v, x * p) * x^v

    y1 = (2 * p)^v * gamma(v + 1/2) / (sqrt(pi) * (z^2 + p^2)^(v + 1/2))
    y2, E = quadde(fun, 0, a)
    Sn = y1 - y2  # "exact" value of the integral tail

    val, E = part_extrap(fun, a, q, z, alpha, tol, kmax, u)
    abs(val - Sn) < tol
```

# References
- MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
"""
function part_extrap(fun::Function, a, q, z, alpha, tol = 1e-6, kmax = 10, u = 1)
    X = zeros(kmax + 2)
    X[1] = a
    partial_sums = zeros(ComplexF64, kmax + 2)
    partial_errors = similar(partial_sums)
    for k = 1 : kmax+1
        X[k + 1] = X[k] + q
        i, e = quadde(fun, X[k], X[k + 1])
        partial_sums[k + 1] = i + partial_sums[k]
        partial_errors[k + 1] = e + partial_errors[k]
    end

    R = partial_sums
    local old0 = 0.0
    local old1 = 0.0
    local val = 0.0
    local error_estimate = 0.0
    for k = 0 : kmax + 1
        if k > 0
            # analytical remainder estimates
            Gk = -exp(-q * z) * (X[k] / X[k + 1])^alpha
        else
            Gk = 0
        end

        # Mosig-Michalski extrapolation
        for j = 1:k
            d = X[k - j + 2] - X[k - j + 1]
            eta = Gk / (1 + u * (j - 1) * d / X[k - j + 1])
            R[k - j + 1] = (R[k - j + 2] - eta * R[k - j + 1]) / (1 - eta)
        end
        val = R[1]

        error_estimate = max(abs(val - old0), abs(val - old1))
        if k > 1 && abs(error_estimate) < tol * abs(val)
            break
        end
        old0 = old1
        old1 = val
    end
    return val, error_estimate
end


# Testing

let N = 20, u = 1
    partial_sum(n::Int) = cumsum([series_term(i) for i = 0:(n - 1)])
    
    series_term(i::Int) = (-1)^i / sqrt(i + 1)
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, 0.604898643421630)

    series_term(i::Int) = (4/5)^(i+1) / (i + 1)
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, log(5))

    series_term(i::Int) = 1 / (i + 1)^2
    extrapolated_sum = mosig_michalski(partial_sum(N), u)
    @assert isapprox(extrapolated_sum, pi^2 / 6; atol = 1e-2)
end


begin
    v = 2
    p = 1
    a = 5.13562
    q = pi / p
    z = 0
    tol = 1e-9
    alpha = 1/2 - v
    kmax = 10
    u = 2

    fun = x -> exp(-x * z) * besselj(v, x * p) * x^v

    y1 = (2 * p)^v * gamma(v + 1/2) / (sqrt(pi) * (z^2 + p^2)^(v + 1/2))
    y2, E = quadde(fun, 0, a)
    Sn = y1 - y2  # "exact" value of the integral tail

    val, E = part_extrap(fun, a, q, z, alpha, tol, kmax, u)
    @assert abs(val - Sn) < tol
end
