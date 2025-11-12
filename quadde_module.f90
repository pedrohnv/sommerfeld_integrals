! Double exponential quadrature integration.
!
! Compile with:
!    gfortran -c quadde_module.f90
!
! # References
! - [1] Bailey, David H., Karthik Jeyabalan, and Xiaoye S. Li. “A comparison of three high-precision quadrature schemes.” Experimental Mathematics 14.3 (2005): 317-329.
! - [2] Vanherck, Joren, Bart Sorée, and Wim Magnus. “Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic.” arXiv preprint arXiv:2007.15057 (2020).
! - [3] van Engelen, Robert A. “Improving the Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas.” https://www.genivia.com/files/qthsh.pdf
! - [4] Takahasi, Hidetosi, and Masatake Mori. “Double exponential formulas for numerical integration.” Publications of the Research Institute for Mathematical Sciences 9.3 (1974): 721-741.
module quadde_module
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: wp => real64  ! working precision
    implicit none

    private  ! All entities are now module-private by default
    public quadde, integrable_function

    abstract interface
        function integrable_function(x) result(y)
            import wp
            real(wp), intent(in) :: x
            real(wp) :: y
        end function integrable_function
    end interface

contains

    !> Double exponential quadrature of f(x) in the interval [a, b] with at most
    !> n levels (more than 6 usually leads to no improvement on convergence)
    !> with requested error eps.
    function quadde(f, a, b, n, eps) result(val)
        procedure(integrable_function) :: f
        integer :: n, k, mode
        real(wp) :: a, b, eps, tol, val, c, d, s, e, h, p, q, fp, fm, v, t, r, x, w, u, eh, y, sign

        if (eps < 0.0) stop "eps must be positive"
        if (n < 0) stop "n must be positive"

        tol = 10.0 * eps
        mode = 0  ! tanh-sinh
        c = 0.0
        d = 1.0
        sign = 1.0
        h = 2.0

        if (b < a) then  ! swap bounds
            v = b
            b = a
            a = v
            sign = -1.0
        end if

        if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
            c = (a + b) * 0.5
            d = (b - a) * 0.5
        else if (ieee_is_finite(a)) then
            mode = 1  ! exp-sinh
            c = a
            v = a + d
        else if (ieee_is_finite(b)) then
            mode = 1  ! exp-sinh
            d = -d
            sign = -sign
            c = b
            v = b + d
        else
            mode = 2  ! sinh-sinh
            v = 0.0
        end if

        s = f(v)
        do k = 0, n
            p = 0.0
            fp = 0.0
            fm = 0.0

            h = h * 0.5
            eh = exp(h)
            t = eh

            if (k > 0) eh = eh * eh

            if (mode == 0) then  ! tanh-sinh
                do
                    u = exp(1.0 / t - t)  ! = exp(-2*sinh(j*h)) = 1/exp(sinh(j*h))^2
                    r = 2.0 * u / (1.0 + u)  ! = 1 - tanh(sinh(j*h))
                    w = (t + 1.0 / t) * r / (1.0 + u)  ! = cosh(j*h) / cosh(sinh(j*h))^2
                    x = d * r

                    if (a + x > a) then  ! too close to a then reuse previous fp
                        y = f(a + x)
                        if (ieee_is_finite(y)) fp = y
                    end if

                    if (b - x < b) then  ! too close to b then reuse previous fm
                        y = f(b - x)
                        if (ieee_is_finite(y)) fm = y
                    end if

                    q = w * (fp + fm)
                    p = p + q
                    t = t * eh

                    if (abs(q) < eps * abs(p)) exit
                end do
            else
                t = t * 0.5
                do
                    r = exp(t - 0.25 / t)  ! exp(sinh(j*h))
                    w = r
                    q = 0.0
                    if (mode == 1) then  ! exp-sinh
                        x = c + d / r
                        if (abs(x - c) < epsilon(x)) exit  ! x hit the finite endpoint
                        y = f(x)
                        if (ieee_is_finite(y)) q = q + y / w
                    else  ! sinh-sinh
                        r = (r - 1.0 / r) * 0.5  ! sinh(sinh(j*h))
                        w = (w + 1.0 / w) * 0.5  ! cosh(sinh(j*h))
                        x = c - d * r
                        y = f(x)
                        if (ieee_is_finite(y)) q = q + y * w
                    end if
                    x = c + d * r
                    y = f(x)
                    if (ieee_is_finite(y)) q = q + y * w
                    q = q * (t + 0.25 / t)  ! cosh(j*h)
                    p = p + q
                    t = t * eh
                    if (abs(q) <= eps * abs(p)) exit
                end do
            end if
            v = s - p
            s = s + p
            if (abs(v) <= tol * abs(s)) exit
        end do
        e = abs(v) / (abs(s) + eps)  ! TODO return or print e
        val = sign * d * s * h  ! result with estimated relative error e
    end function quadde

end module
