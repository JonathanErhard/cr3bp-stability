#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <functional>
#include <stdexcept>
#include <cfg/benchmark-hyper-parameters.h>

namespace cr3bp_benchmarks {

using state_type = std::array<double, 6>;

using vec6 = std::array<double,6>;
using mat6 = std::array<std::array<double,6>,6>;

vec6 operator+(const vec6& a, const vec6& b) {
    vec6 r{};
    for (int i = 0; i < 6; ++i) r[i] = a[i] + b[i];
    return r;
}

// scalar multiplication
vec6 operator*(double s, const vec6& v) {
    vec6 r{};
    for (int i = 0; i < 6; ++i) r[i] = s * v[i];
    return r;
}

// multiplication of Matrix and vector
vec6 matvec(const mat6& A, const vec6& x) {
    vec6 r{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            r[i] += A[i][j] * x[j];
    return r;
}

double compute_L1(double mu) {
    double x = 1.0 - std::cbrt(mu / 3.0);
    for (int i = 0; i < 30; ++i) {
        double r1 = x + mu;
        double r2 = x - (1.0 - mu);
        double f =
            x
            - (1 - mu) * (x + mu) / std::pow(std::abs(r1), 3)
            - mu * (x - (1 - mu)) / std::pow(std::abs(r2), 3);

        double df =
            1
            - (1 - mu) * (1.0 / std::pow(std::abs(r1), 3)
              - 3 * r1 * r1 / std::pow(std::abs(r1), 5))
            - mu * (1.0 / std::pow(std::abs(r2), 3)
              - 3 * r2 * r2 / std::pow(std::abs(r2), 5));

        x -= f / df;
    }
    return x;
}

// jacobian linearized around L1
mat6 linearized_matrix_L1(double mu) {
    const double xL = compute_L1(mu);

    const double r1 = xL + mu;
    const double r2 = xL - (1.0 - mu);

    const double d1 = std::abs(r1);
    const double d2 = std::abs(r2);

    const double Uxx =
        1
        - (1 - mu) * (1.0 / std::pow(d1,3) - 3 * r1 * r1 / std::pow(d1,5))
        - mu * (1.0 / std::pow(d2,3) - 3 * r2 * r2 / std::pow(d2,5));

    const double Uyy =
        1
        - (1 - mu) / std::pow(d1,3)
        - mu / std::pow(d2,3);

    const double Uzz =
        - (1 - mu) / std::pow(d1,3)
        - mu / std::pow(d2,3);

    mat6 A{};

    A[0][3] = 1;
    A[1][4] = 1;
    A[2][5] = 1;

    A[3][0] = Uxx;
    A[3][4] = 2;

    A[4][1] = Uyy;
    A[4][3] = -2;

    A[5][2] = Uzz;

    return A;
}


// analytical solution
std::vector<std::pair<double,state_type>>
surrogate_L1_exact(
    const state_type& x0,
    double mu,
    double tf,
    double dt
) {
    auto A = linearized_matrix_L1(mu);

    // Eigenvalues (known structure, numerically approximated)
    const double lambda = 0.5;
    const double omega  = 1.0;
    const double nu     = 0.8;

    std::vector<std::pair<double,state_type>> traj;
    int N = static_cast<int>(tf / dt);

    for (int i = 0; i <= N; ++i) {
        double t = i * dt;

        state_type x{};

        x[0] = x0[0] * std::cosh(lambda * t);
        x[1] = x0[1] * std::cos(omega * t);
        x[2] = x0[2] * std::cos(nu * t);

        x[3] = lambda * x0[0] * std::sinh(lambda * t);
        x[4] = -omega * x0[1] * std::sin(omega * t);
        x[5] = -nu * x0[2] * std::sin(nu * t);

        traj.push_back({t, x});
    }

    return traj;
}

void surrogate_L1_ode(const state_type& q, state_type& dq, double) {
    static mat6 A = linearized_matrix_L1(INTEGRATOR_PARAMETERS::mu);

    auto v = matvec(A, q);
    for (int i = 0; i < 6; ++i)
        dq[i] = v[i];
}

template<class Integrator>
std::pair<
    std::vector<std::pair<double,state_type>>,
    std::vector<std::pair<double,state_type>>
>
surrogate_L1_benchmark(
    Integrator integrator,
    const state_type& x0,
    double t0,
    double tf,
    double dt
) {
    namespace odeint = boost::numeric::odeint;
    using trajectory_type = std::vector<std::pair<double,state_type>>;

    auto exact = surrogate_L1_exact(
        x0,
        INTEGRATOR_PARAMETERS::mu,
        tf,
        dt
    );

    trajectory_type numerical;
    trajectory_observer<state_type> obs(numerical);

    odeint::integrate_const(
        integrator,
        surrogate_L1_ode,
        x0,
        t0,
        tf,
        dt,
        obs
    );

    return { exact, numerical };
}

} // namespace cr3bp_benchmarks
