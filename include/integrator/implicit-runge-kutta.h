#pragma once
#include <array>
#include <functional>
#include <iostream>
#include <Eigen/Dense>

#include "runge-kutta.h"

/**
 * @brief General implicit Runge-Kutta integrator
 *
 * @tparam T dimension of the state
 * @tparam K number of stages
 */
template <size_t T, size_t K>
class ImplicitRungeKutta : public RungeKutta<T, K> {
public:
// maybe I should use that convetion in the other files as well :/
    using VectorT = Eigen::Matrix<double, T, 1>;
    using MatrixT = Eigen::Matrix<double, T, T>;
    using StageVector = Eigen::Matrix<double, T*K, 1>;
    using StageMatrix = Eigen::Matrix<double, T*K, T*K>;

    ImplicitRungeKutta(
        std::function<VectorT(double, VectorT)> function,
        const Eigen::Matrix<double, K, K>& a,
        const Eigen::Matrix<double, K, 1>& b,
        size_t newton_iters = 10,
        double tol = 1e-10
    )
        : RungeKutta<T, K>(function, a, b),
          m_newton_iters(newton_iters),
          m_tol(tol)
    {}

    void step(double dt) {
        // Stage variables stacked into one vector
        StageVector G = StageVector::Zero();

        // Initial guess: explicit Euler
        for (size_t i = 0; i < K; ++i) {
            G.template segment<T>(i*T) =
                this->m_function(this->m_time, this->m_state);
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
                VectorT fi = this->m_function(ti, stage_sum);

                // Residual: Gi - f(...)
                F.template segment<T>(i*T) =
                    G.template segment<T>(i*T) - fi;

                // Jacobian blocks
                MatrixT df_dx = numerical_jacobian(ti, stage_sum);

                for (size_t j = 0; j < K; ++j) {
                    MatrixT block = MatrixT::Zero();
                    if (i == j) {
                        block.setIdentity();
                    }
                    block -= dt * this->m_A(i, j) * df_dx;

                    J.template block<T,T>(i*T, j*T) = block;
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

    MatrixT numerical_jacobian(double t, const VectorT& x) {
        constexpr double eps = 1e-8;
        MatrixT J;
        VectorT fx = this->m_function(t, x);

        for (size_t i = 0; i < T; ++i) {
            VectorT xp = x;
            xp(i) += eps;
            J.col(i) = (this->m_function(t, xp) - fx) / eps;
        }
        return J;
    }
};
