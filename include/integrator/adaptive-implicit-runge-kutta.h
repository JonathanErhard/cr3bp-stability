#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <Eigen/Dense>

#include "runge-kutta.h"

/*
Adaptive implicit Runge–Kutta integrator.

Step-size adaptation follows the scheme described for adaptive methods and
variable-step implicit RK in Atallah et al. (2020):
  - Use a scaled local error estimate e(h)
  - Choose h* = h * beta * clamp( (tol/e)^(1/(p+1)), minScale, maxScale )
  - For implicit RK, estimate e(h) from the rate of convergence of the internal
    iteration:  e(h) ≈ | y_{n+1}^{(m)} - y_{n+1}^{(m-1)} | / | y_{n+1}^{(m-1)} |

Here we implement the same concept with a robust, tolerance-scaled RMS norm.
We compute y_{n+1}^{(m)} after each Newton iteration and estimate e(h) from the
change between successive iterates. If Newton fails to converge within the allowed
iterations, the final e(h) is used to reduce the step and retry.

Assumes RungeKutta<System,T,K> provides:
  - m_system(t, x) -> VectorT
  - m_time (double), m_state (VectorT)
  - m_A (KxK), m_B (Kx1) and m_C (Kx1) nodes

Notes:
  - This header keeps an embedded low-order weight vector (b_low) for API
    compatibility and optional diagnostics, but the *step-size controller* uses the
    Atallah-style convergence-based error estimate.
  - Jacobian is computed with simple forward differences, as in the original file.
    (If you want the AutoDiff Jacobian version here as well, say so and I will
    merge them.)
*/

template <class System, size_t T, size_t K>
class AdaptiveImplicitRungeKutta : public RungeKutta<System, T, K> {
public:
    using VectorT      = Eigen::Matrix<double, T, 1>;
    using MatrixT      = Eigen::Matrix<double, T, T>;
    using StageVector  = Eigen::Matrix<double, T*static_cast<int>(K), 1>;
    using StageMatrix  = Eigen::Matrix<double, T*static_cast<int>(K), T*static_cast<int>(K)>;

    AdaptiveImplicitRungeKutta(
        System system,
        const Eigen::Matrix<double, K, K>& a,
        const Eigen::Matrix<double, K, 1>& b_high,
        const Eigen::Matrix<double, K, 1>& b_low,
        int order_high,
        size_t newton_iters = 10,
        double atol = 1e-9,
        double rtol = 1e-6,
        double safety_beta = 0.9,
        double min_scale = 0.2,
        double max_scale = 5.0,
        size_t max_rejects = 12
    )
        : RungeKutta<System, T, K>(system, a, b_high),
          m_b_low(b_low),
          m_order_high(order_high),
          m_newton_iters(newton_iters),
          m_atol(atol),
          m_rtol(rtol),
          m_beta(safety_beta),
          m_min_scale(min_scale),
          m_max_scale(max_scale),
          m_max_rejects(max_rejects)
    {}

    // Attempts a step; may internally reject and retry with smaller dt.
    // On success: advances (state,time), updates dt for next call, returns controller error norm (<= 1).
    // On failure: returns NaN (dt is reduced).
    double step(double& dt) {
        if (!(dt > 0.0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        for (size_t attempt = 0; attempt < m_max_rejects; ++attempt) {
            StageVector G = initial_stage_guess();

            double iter_err = std::numeric_limits<double>::infinity();
            const bool newton_ok = solve_stages_newton(G, dt, iter_err);

            // If Newton failed, reject and reduce dt using the same controller formula
            // (Atallah: if convergence not achieved within M, reject and use final e).
            if (!newton_ok) {
                dt *= choose_step_factor(iter_err);
                continue;
            }

            // Accepted-stage solution -> compute high order update (and optional low order one)
            const VectorT y0 = this->m_state;
            const VectorT y_high = y0 + dt * weighted_stage_sum(G, this->m_B);

            // Optional embedded estimate (not used for adaptation here)
            // VectorT y_low = y0 + dt * weighted_stage_sum(G, m_b_low);
            // VectorT embedded_e = y_high - y_low;

            // Controller uses convergence-based estimate (Atallah Eq. 21), scaled with (atol,rtol)
            const double err_norm = std::isfinite(iter_err) ? iter_err : std::numeric_limits<double>::infinity();

            const double factor = choose_step_factor(err_norm);

            if (err_norm <= 1.0) {
                // Accept
                this->m_state = y_high;
                this->m_time += dt;
                dt *= factor;
                return err_norm;
            }

            // Reject
            dt *= factor;
        }

        return std::numeric_limits<double>::quiet_NaN();
    }

    void set_tolerances(double atol, double rtol) { m_atol = atol; m_rtol = rtol; }
    double atol() const { return m_atol; }
    double rtol() const { return m_rtol; }

private:
    Eigen::Matrix<double, K, 1> m_b_low;

    int    m_order_high;
    size_t m_newton_iters;

    double m_atol;
    double m_rtol;

    double m_beta;
    double m_min_scale;
    double m_max_scale;
    size_t m_max_rejects;

    // Scaled RMS error norm (Hairer-style), returns dimensionless value.
    double scaled_rms_norm(const VectorT& err, const VectorT& y0, const VectorT& y1) const {
        VectorT scale;
        for (size_t i = 0; i < T; ++i) {
            const double s = m_atol + m_rtol * std::max(std::abs(y0(static_cast<int>(i))),
                                                       std::abs(y1(static_cast<int>(i))));
            scale(static_cast<int>(i)) = (s > 0.0) ? s : 1.0;
        }
        const VectorT z = err.cwiseQuotient(scale);
        return std::sqrt(z.squaredNorm() / static_cast<double>(T));
    }

    // Atallah-style step factor (Eq. 15), with tol = 1.0 in the scaled norm.
    double choose_step_factor(double err_norm) const {
        const double expo = 1.0 / (static_cast<double>(m_order_high) + 1.0);

        if (!std::isfinite(err_norm)) {
            // catastrophic failure -> reduce aggressively
            return m_min_scale;
        }
        if (err_norm == 0.0) {
            return m_max_scale;
        }

        double s = m_beta * std::pow(1.0 / err_norm, expo);
        return std::clamp(s, m_min_scale, m_max_scale);
    }

    StageVector initial_stage_guess() const {
        StageVector G = StageVector::Zero();
        const VectorT f0 = this->m_system(this->m_time, this->m_state);
        for (size_t i = 0; i < K; ++i) {
            G.template segment<T>(static_cast<int>(i)*static_cast<int>(T)) = f0;
        }
        return G;
    }

    VectorT weighted_stage_sum(const StageVector& G, const Eigen::Matrix<double, K, 1>& weights) const {
        VectorT sum = VectorT::Zero();
        for (size_t i = 0; i < K; ++i) {
            sum += weights(static_cast<int>(i)) *
                   G.template segment<T>(static_cast<int>(i)*static_cast<int>(T));
        }
        return sum;
    }

    // Newton solve for stage derivatives G.
    // Outputs iter_err: convergence-based error estimate using y_{n+1}^{(m)} - y_{n+1}^{(m-1)}.
    bool solve_stages_newton(StageVector& G, double dt, double& iter_err) {
        // Newton convergence tolerance (scaled with dt, atol, rtol)
        const double newton_abs = 0.1 * m_atol / std::max(dt, 1e-12);
        const double newton_rel = 0.1 * m_rtol;

        // Build y_{n+1}^{(0)} from initial guess for convergence-based error estimate
        const VectorT y0 = this->m_state;
        VectorT y_prev = y0 + dt * weighted_stage_sum(G, this->m_B);
        iter_err = std::numeric_limits<double>::infinity();

        for (size_t iter = 0; iter < m_newton_iters; ++iter) {
            StageVector F = StageVector::Zero();
            StageMatrix J = StageMatrix::Zero();

            for (size_t i = 0; i < K; ++i) {
                VectorT xi = this->m_state;
                for (size_t j = 0; j < K; ++j) {
                    xi += dt * this->m_A(static_cast<int>(i), static_cast<int>(j)) *
                          G.template segment<T>(static_cast<int>(j)*static_cast<int>(T));
                }

                const double ti = this->m_time + dt * this->m_C(static_cast<int>(i));
                const VectorT fi = this->m_system(ti, xi);

                // Residual: F_i(G) = G_i - f(t_i, x_i(G))
                F.template segment<T>(static_cast<int>(i)*static_cast<int>(T)) =
                    G.template segment<T>(static_cast<int>(i)*static_cast<int>(T)) - fi;

                const MatrixT df = numerical_jacobian(ti, xi);

                // dF_i/dG_j = I - dt * a_ij * df/dx
                for (size_t j = 0; j < K; ++j) {
                    MatrixT block = MatrixT::Identity();
                    block -= dt * this->m_A(static_cast<int>(i), static_cast<int>(j)) * df;
                    J.template block<T, T>(static_cast<int>(i)*static_cast<int>(T),
                                           static_cast<int>(j)*static_cast<int>(T)) = block;
                }
            }

            Eigen::FullPivLU<StageMatrix> lu(J);
            if (!lu.isInvertible()) {
                // Use last available iter_err to reduce dt
                return false;
            }

            const StageVector delta = lu.solve(-F);
            if (!delta.allFinite()) {
                return false;
            }

            G += delta;

            // Convergence-based error estimate (Atallah Eq. 21 analogue)
            const VectorT y_curr = y0 + dt * weighted_stage_sum(G, this->m_B);
            iter_err = scaled_rms_norm(y_curr - y_prev, y0, y_curr);
            y_prev = y_curr;

            // Newton convergence check on delta
            const double d_inf = delta.cwiseAbs().maxCoeff();
            const double g_inf = G.cwiseAbs().maxCoeff();
            if (d_inf <= newton_abs + newton_rel * g_inf) {
                // Make sure residual is finite (avoid accepting a diverged residual)
                const double f_inf = F.cwiseAbs().maxCoeff();
                return std::isfinite(f_inf);
            }
        }

        // Did not converge in allotted iterations; iter_err holds the last change in y_{n+1}
        return false;
    }

    MatrixT numerical_jacobian(double t, const VectorT& x) {
        constexpr double eps = 1e-8;
        MatrixT J;
        const VectorT fx = this->m_system(t, x);

        for (size_t i = 0; i < T; ++i) {
            VectorT xp = x;
            xp(static_cast<int>(i)) += eps;
            J.col(static_cast<int>(i)) = (this->m_system(t, xp) - fx) / eps;
        }
        return J;
    }
};