// Variable-step Gauss–Legendre (VGL) integrator for the CR3BP
//  this method is only symplectic for fixed step sizes :/ idk i'll still use it I guess
//
// References:
//  Hairer, Lubich, Wanner — Geometric Numerical Integration
//  Butcher — Numerical Methods for ODEs

// Variable-step Gauss–Legendre (VGL) integrator for the CR3BP 

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

using state_type = Eigen::Matrix<double,6,1>;

// One fixed-step Gauss–Legendre (order 4) step
void gl4_step(const CR3BP &sys, state_type &x, double h)
{
    const double a11 = 0.25;
    const double a12 = 0.25 - std::sqrt(3.0)/6.0;
    const double a21 = 0.25 + std::sqrt(3.0)/6.0;
    const double a22 = 0.25;
    const double b1  = 0.5;
    const double b2  = 0.5;

    state_type k1 = sys(x);
    state_type k2 = k1;

    for (int it = 0; it < 8; ++it) {
        state_type y1 = x + h * (a11 * k1 + a12 * k2);
        state_type y2 = x + h * (a21 * k1 + a22 * k2);
        k1 = sys(y1);
        k2 = sys(y2);
    }

    x += h * (b1 * k1 + b2 * k2);
}

// Adaptive variable-step Gauss–Legendre using step doubling
void adaptive_gl4_step(const CR3BP &sys,
                        state_type &x,
                        double &h,
                        double tol)
{
    const double safety = 0.9;
    const double p = 4.0;

    while (true) {
        state_type x_full = x;
        gl4_step(sys, x_full, h);

        state_type x_half = x;
        gl4_step(sys, x_half, 0.5 * h);
        gl4_step(sys, x_half, 0.5 * h);

        double err = (x_half - x_full).norm();

        if (err < tol || h < 1e-12) {
            x = x_half;
            if (err > 1e-15)
                h *= safety * std::pow(tol / err, 1.0 / (p + 1.0));
            break;
        } else {
            h *= safety * std::pow(tol / err, 1.0 / (p + 1.0));
        }
    }
}

int main()
{
    CR3BP sys{0.0121505856};

    state_type x;
    x << 1.0 - sys.mu - 0.01, 0, 0,
         0, 0.1, 0;

    double t = 0.0;
    double tf = 10.0;
    double h = 0.01;
    double tol = 1e-10;

    while (t < tf) {
        adaptive_gl4_step(sys, x, h, tol);
        t += h;
        std::cout << t << " " << x.transpose() << '\n';
    }
}
