#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>

#include <eigen3/Eigen/Dense>

#include "integrator/adaptive-implicit-runge-kutta.h"
#include "model/CR3BP.h"

#include <eigen3/Eigen/Dense>

int main()
{
    constexpr size_t T = 6;
    constexpr size_t K = 2;

    using Vector6 = Eigen::Matrix<double,6,1>;

    const double mu = 0.0121505856;

    // ---- Gauss–Legendre 2-stage (order 4) ----
    Eigen::Matrix<double,K,K> A;
    Eigen::Matrix<double,K,1> b_high, b_low;

    const double s3 = std::sqrt(3.0);

    A << 0.25,             0.25 - s3/6.0,
         0.25 + s3/6.0,    0.25;

    b_high << 0.5, 0.5;

    // simple embedded lower-order rule (not symplectic, but OK for testing)
    b_low << 1.0, 0.0;

    // ---- ODE wrapper ----
    auto rhs = [mu](double /*t*/, Vector6 x) -> Vector6 {
        Vector6 dx;
        CR3BP::calculate_derivative_rotating(x, mu, dx);
        return dx;
    };

    AdaptiveImplicitRungeKutta<T,K> integrator(
        rhs,
        A,
        b_high,
        b_low,
        8,        // Newton iterations
        1e-10,    // tolerance
        0.9       // safety
    );

    // ---- Initial condition (near L4-type test state) ----
    Vector6 state;
    state << 0.48784941, 0.86602540, 0.0,
             -0.86602540, 0.48784941, 0.0;

    integrator.set_state(state);
    integrator.set_time(0.0);

    double dt = 1e-2;
    const double t_end = 5.0;

    Vector6 state_rot = state;
    double H0 = CR3BP::calculate_hamiltonian_rotating(state_rot, mu);

    std::cout << std::setprecision(10);

    while (integrator.time() < t_end) {
        integrator.step(dt);

        Vector6 x = integrator.state();
        double t  = integrator.time();

        Vector6 xrot = x;
        double H = CR3BP::calculate_hamiltonian_rotating(xrot, mu);

        std::cout
            << t << " "
            << x.transpose() << " "
            << (H - H0) << "\n";
    }

    return 0;
}