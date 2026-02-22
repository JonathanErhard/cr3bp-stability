#include <array>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/AutoDiff>
#include "model/CR3BP.h"

namespace cr3bp_benchmarks {
using state_type = std::array<double, 6>;
using Vector6d = Eigen::Matrix<double, 6, 1>;

// Helper methods to convert std::array to Eigen and back
inline Vector6d to_eigen(const state_type& s) { return Eigen::Map<const Vector6d>(s.data()); }
inline void to_array(const Vector6d& v, state_type& s) { Eigen::Map<Vector6d>(s.data()) = v; }

// Matrix exponential using Eigen's ComplexEigenSolver
// I have used a textbook from my first semester as a reference I am guessing that this approach can be found in any linear algebra book
template <typename SCALAR_TYPE>
Eigen::Matrix<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic> matrixExp(const Eigen::Matrix<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic>& A, double t)
{
    Eigen::ComplexEigenSolver<Eigen::Matrix<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic>> ces(A);
    Eigen::Matrix<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic> V = ces.eigenvectors();
    Eigen::Matrix<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic> D = ces.eigenvalues().asDiagonal();

    // Exponentiate eigenvalues 
    for (int i = 0; i < D.rows(); i++)
        D(i,i) = exp(D(i,i) * t);

    // Reconstruct e^(At)
    Eigen::MatrixXcd expAt = V * D * V.inverse();

    return expAt.real();
}

// Jacobian calculation using Eigen AutoDiff
Eigen::Matrix<double,6,6> calculate_linearized_A(const Vector6d& state, double mu) {
    using ADScalar = Eigen::AutoDiffScalar<Vector6d>;
    using Vector6AD = Eigen::Matrix<ADScalar, 6, 1>;

    Vector6AD state_ad;
    for (int i = 0; i < 6; ++i) {
        state_ad(i).value() = state(i);
        state_ad(i).derivatives() = Eigen::Matrix<double,1,6>::Unit(i);
    }

    Vector6AD deriv_ad;
    CR3BP::calculate_derivative_rotating<ADScalar>(state_ad, mu, deriv_ad);

    Eigen::Matrix<double,6,6> A;
    for (int i = 0; i < 6; ++i) {
        A.row(i) = deriv_ad(i).derivatives();
    }
    return A;
}

// Compute L1 position using Newton's method
double compute_L1(double mu) {
    double x = 1.0 - std::cbrt(mu / 3.0);
    for (int i = 0; i < 10; ++i) {
        Vector6d s = Vector6d::Zero();
        s[0] = x;
        Vector6d ds;
        CR3BP::calculate_derivative_rotating(s, mu, ds);
        
        // Use AutoDiff for the 1D Newton step specifically for L1
        auto A = calculate_linearized_A(s, mu);
        x -= ds[3] / A(3, 0); 
    }
    return x;
}

// Surrogate L1 dynamics ODE to use for numeric integration
void surrogate_L1_ode(const state_type& q, state_type& dq, double) {
    const double mu = INTEGRATOR_PARAMETERS::mu;
    static double xL = compute_L1(mu);
    static Eigen::Matrix<double,6,6> A = calculate_linearized_A((Vector6d() << xL, 0, 0, 0, 0, 0).finished(), mu);

    Vector6d res = A * to_eigen(q);
    to_array(res, dq);
}

// Generate exact surrogate L1 trajectory for benchmarking as a reference trajectory
std::vector<std::pair<double, state_type>>
surrogate_L1_exact(const state_type& x0, double mu, double tf, double dt) {
    const double lambda = 0.5; // Example value
    const double omega = 1.0;
    const double nu = 0.8;

    std::vector<std::pair<double, state_type>> traj;
    for (double t = 0; t <= tf; t += dt) {
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
}