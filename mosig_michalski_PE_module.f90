!> Uso do algoritmo de Mosig-Michalski para computar integrais de Sommerfeld.
!>
!> Dada a seguinte integral que considera uma corrente em uma região circular de
!> raio $\tau$, onde $z$ e $z'$ são distâncias perpendiculares à estratificação,
!> do campo e da fonte, respectivamente, e $\rho$ é a distância lateral entre esses
!> pontos.
!>
!> # Math: \int_0^\infty \tilde{G}(k_\rho; z|z') \, J_\nu(k_\rho \rho) \, J_\sigma(k_\rho \tau) \, k_\rho \, \mathrm{d}k_\rho
!>
!> Estas integrais são frequentemente oscilatórias e divergentes com polos
!> (singularidades) no eixo real em $a_1$, $a_2$, etc. Para computar a cauda desta
!> integral mais rapidamente, faz-se partição e extrapolação (via algoritmo de
!> Mosig-Michalski) da série infinita que aproxima a integral.
!>
!> # Referências
!> - MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
!> - MICHALSKI, K. A. On the efficient evaluation of integral arising in the sommerfeld halfspace problem. In: IEE Proceedings H (Microwaves, Antennas and Propagation). IET Digital Library, 1985. p. 312-318.
module Mosig_Michalski_PE_module
    use iso_fortran_env, only: wp => real64  ! working precision
    use quadde_module  ! quadde, integrable_function
    implicit none

contains

    !> Sequence acceleation of an infinite series by extrapolation via the Mosig–Michalski algorithm (t-transformation weighted averages).
    !>
    !> # Arguments
    !> - partial_sums: Vector of partial sums of the sequence to be accelerated, where the n-th entry should be the partial sum up to index n.
    !> - u: Transformation parameter controlling the extrapolation's sensitivity (1 for logarithmic convergence, 2 otherwise).
    !>
    !> # Returns
    !> - extrapolated_sum: The extrapolated sum, providing an improved estimate for the series true value beyond available terms.
    !>
    !> # Examples
    !> ```Fortran
    !> use Mosig_Michalski_PE
    !> implicit none
    !> real(kind=wp), parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286198_wp
    !> integer, parameter :: N = 20
    !> real(kind=wp) :: u
    !> complex(kind=wp) :: series_term, partial_sums(N), extrapolated_sum
    !> integer :: i
    !>
    !> u = 1.0_wp
    !> extrapolated_sum = 0.0_wp
    !>
    !> do i = 0, N-1
    !>     series_term = (-1.0_wp)**i / sqrt(real(i+1))
    !>     if (i == 0) then
    !>         partial_sums(i + 1) = series_term
    !>     else
    !>         partial_sums(i + 1) = partial_sums(i) + series_term
    !>     end if
    !> end do
    !> extrapolated_sum = mosig_michalski(partial_sums, u)
    !> ```
    !>
    !> # References
    !> - MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
    function mosig_michalski(partial_sums, u) result(extrapolated_sum)
        complex(kind=wp), dimension(:), intent(in) :: partial_sums
        real(kind=wp), intent(in) :: u
        real(kind=wp) :: beta
        integer :: kmax, k, j
        complex(kind=wp), allocatable :: R(:), X(:)
        complex(kind=wp) :: wk0, wk1, Gk, d, eta, extrapolated_sum

        kmax = size(partial_sums) - 1
        allocate(R(kmax + 1))
        allocate(X(kmax + 1))

        beta = 1.0_wp
        do k = 0, kmax
            X(k+1) = k + beta  ! equidistant nodes
            R(k+1) = partial_sums(k+1)
        end do

        ! Transformed sequence values at each extrapolation order k
        ! R = partial_sums

        ! remainder estimates
        wk0 = 1.0_wp  ! w_{k-1}
        wk1 = 1.0_wp  ! w_{k}

        do k = 0, kmax
            if (k >= 1) then
                wk1 = partial_sums(k + 1) - partial_sums(k)
                Gk = wk1 / wk0
            else
                Gk = 0.0_wp
            end if
            do j = 1, k
                d = X(k - j + 2) - X(k - j + 1)
                eta = Gk / (1.0_wp + u * (j - 1) * d / X(k - j + 1))
                R(k - j + 1) = (R(k - j + 2) - eta * R(k - j + 1)) / (1.0_wp - eta)
            end do
            extrapolated_sum = R(1)
            wk0 = wk1
        end do

        deallocate(R)
        deallocate(X)
    end function


    !> Numerical integragion of Sommerfeld integral tails with Partition-Extrapolation
    !> via the Mosig-Michalski method.
    !>
    !> The expected form of the integral is the following, where ``\\tau`` is the
    !> radius of the circular region that contains nonzero current, ``z`` and ``z'``
    !> are the perpendicular distances to the medium stratification of the field
    !> and the source, respectively, and ``\\rho`` is the lateral distance between
    !> those points.
    !>
    !> ```math
    !> \\int_0^\\infty \\tilde{G}(k_\\rho; z|z') \\, J_\\nu(k_\\rho \\rho) \\, J_\\sigma(k_\\rho \\tau) \\, k_\\rho \\, \\mathrm{d}k_\\rho
    !> ```
    !>
    !> # Parameters
    !> - f: The integrand function ``f(x)`` to be integrated over each partition subinterval.
    !> - a: The lower limit of integration (typically the start of the tail, after any head subtraction).
    !> - q: Partition step size; the width of each subinterval (typically matched to Bessel oscillation period or chosen for balance between cost and convergence).
    !> - z: Tail asymptotic decay parameter (often related to physical geometry or stratification; used in remainder estimate).
    !> - alpha: Decay exponent in analytic remainder estimate for the extrapolation transformation.
    !> - tol: Desired relative tolerance for extrapolation stopping criterion.
    !> - kmax: Maximum number of partitions/subintervals for extrapolation.
    !> - u: Transformation parameter controlling the extrapolation's sensitivity (1 for logarithmic convergence, 2 otherwise).
    !>
    !> # Returns
    !> - val: Accelerated estimate of the semi-infinite tail integral.
    !> - error_estimate: Estimate of extrapolation absolute error.
    !>
    !> # References
    !> - MICHALSKI, Krzysztof A.; MOSIG, Juan R. Efficient computation of Sommerfeld integral tails – methods and algorithms. Journal of Electromagnetic Waves and Applications, v. 30, n. 3, p. 281-317, 2016.
    subroutine part_extrap(f, a, q, z, alpha, tol, kmax, u, val, error_estimate)
        procedure(integrable_function) :: f
        real(kind=wp), intent(in) :: a, q, z, alpha, tol, u
        integer, intent(in) :: kmax
        complex(kind=wp), intent(out) :: val, error_estimate
        real(kind=wp) :: d
        integer :: N, k, j
        real(kind=wp), allocatable :: X(:)
        complex(kind=wp), allocatable :: R(:), partial_sums(:)
        complex(kind=wp) :: old0, old1, Gk, eta

        N = kmax + 2
        allocate(X(N))
        allocate(partial_sums(N))
        X(1) = a
        partial_sums(1) = (0.0, 0.0)

        do k = 1, N-1
            X(k + 1) = X(k) + q
            val = quadde(f, X(k), X(k + 1), 6, tol)
            partial_sums(k + 1) = val + partial_sums(k)
        end do

        R = partial_sums
        old0 = (0.0, 0.0)
        old1 = (0.0, 0.0)
        val = (0.0, 0.0)
        error_estimate = (0.0, 0.0)

        do k = 0, kmax + 1
            if (k > 0) then
                ! analytical remainder estimates
                Gk = -exp(-q * z) * (X(k) / X(k + 1))**alpha
            else
                Gk = 0.0
            end if

            ! Mosig-Michalski extrapolation
            do j = 1, k
                d = X(k - j + 2) - X(k - j + 1)
                eta = Gk / (1.0 + u * (j - 1.0) * d / X(k - j + 1))
                R(k - j + 1) = (R(k - j + 2) - eta * R(k - j + 1)) / (1.0 - eta)
            end do
            val = R(1)

            error_estimate = max(abs(val - old0), abs(val - old1))
            if (k > 1 .and. abs(error_estimate) < tol * abs(val)) then
                exit
            end if
            old0 = old1
            old1 = val
        end do

        deallocate(X)
        deallocate(R)
        deallocate(partial_sums)
    end subroutine

end module Mosig_Michalski_PE_module
