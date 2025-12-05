// Compile and run with `gcc -o a.out quadde.c -lm && ./a.out`
//
// # References:
// van Engelen, Robert A. “Improving the Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas.” https://www.genivia.com/files/qthsh.pdf


#include "stdio.h"
#include "float.h"
#include "math.h"

// integrate function f, range a..b, max levels n, error tolerance eps
double quadde(double (*f)(double), double a, double b, int n, double eps) {
    const double tol = 10*eps;
    double c = 0, d = 1, s, sign = 1, e, v, h = 2;
    int k = 0, mode = 0; // Tanh-Sinh = 0, Exp-Sinh = 1, Sinh-Sinh = 2
    if (b < a) { // swap bounds
        v = b;
        b = a;
        a = v;
        sign = -1;
    }
    if (isfinite(a) && isfinite(b)) {
        c = (a+b)/2;
        d = (b-a)/2;
        v = c;
    }
    else if (isfinite(a)) {
        mode = 1; // Exp-Sinh
        // alternatively d = exp_sinh_opt_d(f, a, eps, d);
        c = a;
        v = a+d;
    }
    else if (isfinite(b)) {
        mode = 1; // Exp-Sinh
        d = -d; // alternatively d = exp_sinh_opt_d(f, b, eps, -d);
        sign = -sign;
        c = b;
        v = b+d;
    }
    else {
        mode = 2; // Sinh-Sinh
        v = 0;
    }
    s = f(v);
    do {
        double p = 0, q, fp = 0, fm = 0, t, eh;
        h /= 2;
        t = eh = exp(h);
        if (k > 0)
            eh *= eh;
        if (mode == 0) { // Tanh-Sinh
            do {
                double u = exp(1/t-t); // = exp(-2*sinh(j*h)) = 1/exp(sinh(j*h))^2
                double r = 2*u/(1+u); // = 1 - tanh(sinh(j*h))
                double w = (t+1/t)*r/(1+u); // = cosh(j*h)/cosh(sinh(j*h))^2
                double x = d*r;
                if (a+x > a) { // if too close to a then reuse previous fp
                    double y = f(a+x);
                    if (isfinite(y))
                        fp = y; // if f(x) is finite, add to local sum
                }
                if (b-x < b) { // if too close to a then reuse previous fp
                    double y = f(b-x);
                    if (isfinite(y))
                        fm = y; // if f(x) is finite, add to local sum
                }
                q = w*(fp+fm);
                p += q;
                t *= eh;
            } while (fabs(q) > eps*fabs(p));
        }
        else {
            t /= 2;
            do {
                double r = exp(t-.25/t); // = exp(sinh(j*h))
                double x, y, w = r;
                q = 0;
                if (mode == 1) { // Exp-Sinh
                    x = c + d/r;
                    if (fabs(x - c) < DBL_EPSILON) // if x hit the finite endpoint then break
                        break;
                    y = f(x);
                    if (isfinite(y)) // if f(x) is finite, add to local sum
                        q += y/w;
                }
                else { // Sinh-Sinh
                    r = (r-1/r)/2; // = sinh(sinh(j*h))
                    w = (w+1/w)/2; // = cosh(sinh(j*h))
                    x = c - d*r;
                    y = f(x);
                    if (isfinite(y)) // if f(x) is finite, add to local sum
                        q += y*w;
                }
                x = c + d*r;
                y = f(x);
                if (isfinite(y)) // if f(x) is finite, add to local sum
                    q += y*w;
                q *= t+.25/t; // q *= cosh(j*h)
                p += q;
                t *= eh;
            } while (fabs(q) > eps*fabs(p));
        }
        v = s-p;
        s += p;
        ++k;
    } while (fabs(v) > tol*fabs(s) && k <= n);
    e = fabs(v)/(fabs(s)+eps);
    return sign*d*s*h; // result with estimated relative error e
}



double f1(double x) {
    return 1.0 / ( (x - 2.0) * pow(1.0 - x, 0.25) * pow(1.0 + x, 0.75) );
}

double f2(double x) {
    return cos(M_PI * x) / sqrt(1.0 - x);
}

double f3(double x) {
    return exp(-1.0 - x) / (1.0 + x);
}

double f4(double x) {
    return 1.0 / pow(1.0 + pow(x, 2.0), (5.0 / 4.0));
}

double f5(double x) {
    return 1.0 / (1.0 + pow(x, 4.0));
}

double f6(double x) {
    return 1.0 / pow(x, 2.0);
}


int main() {
    double tol = 1e-6;
    double val;

    val = quadde(f1, -1.0, 1.0, 6, tol);
    printf("Test 1: %.17g\n", val);
    if (fabs(-1.9490 - val) > 1e-4) {
        printf("Test 1 failed");
        return 1;
    }

    val = quadde(f2, -1.0, 1.0, 6, tol);
    printf("Test 2: %.17g\n", val);
    if (fabs(-0.69049 - val) > 1e-4) {
        return 2;
    }

    val = quadde(f3, 0.0, INFINITY, 6, tol);
    printf("Test 3: %.17g\n", val);
    if (fabs(0.21938 - val) > 1e-4) {
        return 3;
    }

    val = quadde(f4, -INFINITY, INFINITY, 6, tol);
    printf("Test 4: %.17g\n", val);
    if (fabs(2.3962 - val) > 1e-4) {
        return 4;
    }

    val = quadde(f5, -INFINITY, INFINITY, 6, tol);
    printf("Test 5: %.17g\n", val);
    if (fabs(2.2214 - val) > 1e-4) {
        return 5;
    }
    
    val = quadde(f6, 0.1, 1.0, 6, tol);
    printf("Test 6: %.17g\n", val);
    if (fabs(9.0 - val) > 1e-4) {
        return 6;
    }

    printf("All Tests passed!\n");

    return 0;
}
