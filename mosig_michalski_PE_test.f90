! Test Mosig-Michalski algorithm for computing Sommerfeld integrals.
program mosig_michalski_PE_test
    use, intrinsic :: ieee_arithmetic
    use Mosig_Michalski_PE_module
    implicit none

    real(kind=wp), parameter :: PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286198_wp
    integer, parameter :: N = 20
    real(kind=wp) :: v, p, a, q, tol, alpha, u, z
    complex(kind=wp) :: series_term, partial_sums(N), extrapolated_sum, val, error_estimate
    integer :: i, kmax
    logical :: is_nan

    print *, "Testing mosig_michalski_PE..."

    ! Series acceleration
    u = 1.0
    extrapolated_sum = (0.0, 0.0)

    ! Example 1
    do i = 0, N-1
        series_term = (-1.0_wp)**real(i) / sqrt(real(i+1))
        if (i == 0) then
            partial_sums(i + 1) = series_term
        else
            partial_sums(i + 1) = partial_sums(i) + series_term
        end if
    end do
    extrapolated_sum = mosig_michalski(partial_sums, u)
    is_nan = ieee_is_nan(real(extrapolated_sum)) .or. ieee_is_nan(aimag(extrapolated_sum))
    if (abs(real(extrapolated_sum) - 0.604898643421630_wp) > 1e-6 .or. is_nan) stop "Test 1 failed"

    ! Example 2
    do i = 0, N-1
        series_term = (4.0_wp / 5.0_wp)**real(i+1) / real(i+1)
        if (i == 0) then
            partial_sums(i + 1) = series_term
        else
            partial_sums(i + 1) = partial_sums(i) + series_term
        end if
    end do
    extrapolated_sum = mosig_michalski(partial_sums, u)
    is_nan = ieee_is_nan(real(extrapolated_sum)) .or. ieee_is_nan(aimag(extrapolated_sum))
    if (abs(real(extrapolated_sum) - log(5.0_wp)) > 1e-6 .or. is_nan) stop "Test 2 failed"

    ! Example 3
    do i = 0, N-1
        series_term = 1.0_wp / real(i+1)**2.0_wp
        if (i == 0) then
            partial_sums(i + 1) = series_term
        else
            partial_sums(i + 1) = partial_sums(i) + series_term
        end if
    end do
    extrapolated_sum = mosig_michalski(partial_sums, u)
    is_nan = ieee_is_nan(real(extrapolated_sum)) .or. ieee_is_nan(aimag(extrapolated_sum))
    if (abs(real(extrapolated_sum) - (PI**2.0_wp / 6.0_wp)) > 1e-2 .or. is_nan) stop "Test 3 failed"

    ! Partition-Extrapolation
    v = 2.0
    p = 1.0
    a = 5.13562
    q = PI / p
    z = 0.0
    tol = 1e-9
    alpha = 1.0 / 2.0 - v
    kmax = 10
    u = 2.0
    ! It is safe to call f(x) onward
    call part_extrap(f, a, q, z, alpha, tol, kmax, u, val, error_estimate)
    print *, "I =", val
    print *, "E =", error_estimate
    is_nan = ieee_is_nan(real(val)) .or. ieee_is_nan(aimag(val))
    if (abs(-10.07948621951323 - val) > 1e-6 .or. is_nan) stop "Partition-Extrapolation failed"

    print *, "All tests passed!"

contains

    function f(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: ZR, ZI, CYR, CYI, FNU
        complex(wp) :: y, besselj
        integer :: IERR, NT, NZ, KODE
        NT = 1
        NZ = 1
        KODE = 1  ! no scaling
        FNU = v
        ZR = x
        ZI = 0.0
        IERR = 0
        call ZBESJ(ZR, ZI, FNU, KODE, NT, CYR, CYI, NZ, IERR)
        if (IERR /= 0) then
            print *, "IERR =", IERR
        end if
        besselj = cmplx(CYR, CYI, kind=wp)
        y = exp(-x * z) * besselj * x**v
    end function

end program mosig_michalski_PE_test
