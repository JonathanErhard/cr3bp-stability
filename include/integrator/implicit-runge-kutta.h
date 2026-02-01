#pragma once
#include <array>
#include <functional>
#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

#include "runge-kutta.h"

/**
 * @brief General implicit Runge-Kutta integrator (Newton solve on stages)
 *
 * Expectations on `System`:
 *   - `System` must be callable with both double-valued and AutoDiff-valued states.
 *   - In particular, `system(t, x)` should work when `x` is:
 *       Eigen::Matrix<double, T, 1> and Eigen::Matrix<ADScalar, T, 1>
 *
 * This allows computing the stage Jacobian via automatic differentiation
 * (no finite differences).
 *
 * @tparam System RHS functor/system (templated/overloaded to accept AutoDiff scalar)
 * @tparam T      dimension of the state
 * @tparam K      number of stages
 */
template <class System, size_t T, size_t K>
class ImplicitRungeKutta : public RungeKutta<System, T, K> {
public:
    using VectorT      = Eigen::Matrix<double, T, 1>;
    using MatrixT      = Eigen::Matrix<double, T, T>;
    using StageVector  = Eigen::Matrix<double, T*K, 1>;
    using StageMatrix  = Eigen::Matrix<double, T*K, T*K>;

    // AutoDiff types for Jacobian calculation
    using DerivVector  = Eigen::Matrix<double, T, 1>;
    using ADScalar     = Eigen::AutoDiffScalar<DerivVector>;
    using VectorAD     = Eigen::Matrix<ADScalar, T, 1>;

    ImplicitRungeKutta(
        System system,
        const Eigen::Matrix<double, K, K>& a,
        const Eigen::Matrix<double, K, 1>& b,
        size_t newton_iters = 10,
        double tol = 1e-10
    )
        : RungeKutta<System, T, K>(system, a, b),
          m_newton_iters(newton_iters),
          m_tol(tol)
    {}

    void step(double dt) {
        // Stage variables stacked into one vector
        StageVector G = StageVector::Zero();

        // Initial guess: explicit Euler on the stage derivatives
        for (size_t i = 0; i < K; ++i) {
            G.template segment<T>(i*T) =
                this->m_system(this->m_time, this->m_state);
        }

        // Newton iteration
        for (size_t iter = 0; iter < m_newton_iters; ++iter) {
            StageVector F = StageVector::Zero();
            StageMatrix J = StageMatrix::Zero();

            for (size_t i = 0; i < K; ++i) {
                VectorT stage_sum = this->m_state;

                for (size_t j = 0; j < K; ++j) {
                    stage_sum += dt * this->m_A(i, j)
                               * G.template segment<T>(j*T);
                }

                double ti = this->m_time + dt * this->m_C(i);
                VectorT fi = this->m_system(ti, stage_sum);

                // Residual: Gi - f(ti, stage_sum)
                F.template segment<T>(i*T) =
                    G.template segment<T>(i*T) - fi;

                // Jacobian df/dx at (ti, stage_sum) via AutoDiff
                MatrixT df_dx = autodiff_jacobian(ti, stage_sum);

                // Assemble block Jacobian for the stage system
                for (size_t j = 0; j < K; ++j) {
                    MatrixT block = MatrixT::Zero();
                    if (i == j) {
                        block.setIdentity();
                    }
                    block -= dt * this->m_A(i, j) * df_dx;

                    J.template block<T, T>(i*T, j*T) = block;
                }
            }

            // Newton update
            StageVector delta = J.fullPivLu().solve(-F);
            G += delta;

            if (delta.norm() < m_tol) {
                break;
            }
        }

        // Final state update
        VectorT derivative = VectorT::Zero();
        for (size_t i = 0; i < K; ++i) {
            derivative += this->m_B(i)
                        * G.template segment<T>(i*T);
        }

        this->m_state += dt * derivative;
        this->m_time  += dt;
    }

private:
    size_t m_newton_iters;
    double m_tol;

    MatrixT autodiff_jacobian(double t, const VectorT& x) {
        // Seed AutoDiff variables with identity derivatives
        VectorAD x_ad;
        for (size_t i = 0; i < T; ++i) {
            x_ad(i).value() = x(static_cast<int>(i));
            x_ad(i).derivatives() = DerivVector::Zero();
            x_ad(i).derivatives()(static_cast<int>(i)) = 1.0;
        }

        // Evaluate system in AutoDiff mode
        VectorAD f_ad = this->m_system(t, x_ad);

        // Extract Jacobian: J(r,c) = d f_r / d x_c
        MatrixT J;
        for (size_t r = 0; r < T; ++r) {
            J.row(static_cast<int>(r)) =
                f_ad(static_cast<int>(r)).derivatives().transpose();
        }
        return J;
    }
};