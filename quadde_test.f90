! Double exponential quadrature integration Test.
!
! # References
! - [1] Bailey, David H., Karthik Jeyabalan, and Xiaoye S. Li. “A comparison of three high-precision quadrature schemes.” Experimental Mathematics 14.3 (2005): 317-329.
! - [2] Vanherck, Joren, Bart Sorée, and Wim Magnus. “Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic.” arXiv preprint arXiv:2007.15057 (2020).
! - [3] van Engelen, Robert A. “Improving the Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas.” https://www.genivia.com/files/qthsh.pdf
! - [4] Takahasi, Hidetosi, and Masatake Mori. “Double exponential formulas for numerical integration.” Publications of the Research Institute for Mathematical Sciences 9.3 (1974): 721-741.

program test_quadde
    use quadde_module
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: wp => real64  ! working precision
    implicit none

    real(wp), parameter :: tol = 1e-6
    real(wp), parameter :: PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286198_wp
    real(wp) :: val, POS_INF, NEG_INF

    POS_INF = ieee_value(1.0_wp, ieee_positive_inf)
    NEG_INF = ieee_value(-1.0_wp, ieee_negative_inf)

    print *, "Testing quadde..."

    val = quadde(f1, -1.0_wp, 1.0_wp, 6, tol)
    print *, "Test 1:", val
    if (abs(-1.9490 - val) > 1e-4 .or. ieee_is_nan(val)) stop "Test 1 failed"

    val = quadde(f2, -1.0_wp, 1.0_wp, 6, tol)
    print *, "Test 2:", val
    if (abs(-0.69049 - val) > 1e-4 .or. ieee_is_nan(val)) stop "Test 2 failed"

    val = quadde(f3, 0.0_wp, POS_INF, 6, tol)
    print *, "Test 3:", val
    if (abs(0.21938 - val) > 1e-4) stop "Test 3 failed"

    val = quadde(f4, NEG_INF, POS_INF, 6, tol)
    print *, "Test 4:", val
    if (abs(2.3962 - val) > 1e-4 .or. ieee_is_nan(val)) stop "Test 4 failed"

    val = quadde(f5, NEG_INF, POS_INF, 6, tol)
    print *, "Test 5:", val
    if (abs(2.2214 - val) > 1e-4 .or. ieee_is_nan(val)) stop "Test 5 failed"

    print *, "All tests passed!"

contains
    ! After contains, you can only define subroutines or functions.
    ! You cannot put assignments, if statements, or any "straight-line" code there.

    function f1(x) result(y)
        real(wp), intent(in) :: x
        complex(wp) :: y
        y = 1.0 / ( (x - 2.0) * (1.0 - x)**0.25 * (1.0 + x)**0.75 )
    end function

    
    function f2(x) result(y)
        real(wp), intent(in) :: x
        complex(wp) :: y
        y = cos(PI * x) / sqrt(1.0 - x)
    end function


    function f3(x) result(y)
        real(wp), intent(in) :: x
        complex(wp) :: y
        y = exp(-1.0 - x) / (1.0 + x)
    end function


    function f4(x) result(y)
        real(wp), intent(in) :: x
        complex(wp) :: y
        y = 1.0 / (1.0 + x**2.0)**(5.0 / 4.0)
    end function

    
    function f5(x) result(y)
        real(wp), intent(in) :: x
        complex(wp) :: y
        y = 1.0 / (1.0 + x**4.0)
    end function

end program test_quadde
