#pragma once
#include <array>
#include <Eigen/Dense>

#include "runge-kutta.h"

template <size_t T, size_t K>
class AdaptiveImplicitRungeKutta : public RungeKutta<T, K> {
public:
    using VectorT = Eigen::Matrix<double, T, 1>;
    using MatrixT = Eigen::Matrix<double, T, T>;
    using StageVector = Eigen::Matrix<double, T*K, 1>;
    using StageMatrix = Eigen::Matrix<double, T*K, T*K>;

    AdaptiveImplicitRungeKutta(
        std::function<VectorT(double, VectorT)> function,
        const Eigen::Matrix<double, K, K>& a,
        const Eigen::Matrix<double, K, 1>& b_high,
        const Eigen::Matrix<double, K, 1>& b_low,
        size_t newton_iters = 8,
        double tol = 1e-8,
        double safety = 0.9
    )
        : RungeKutta<T, K>(function, a, b_high),
          m_b_low(b_low),
          m_newton_iters(newton_iters),
          m_tol(tol),
          m_safety(safety)
    {}

    double step(double& dt) {
        StageVector G = StageVector::Zero();

        for (size_t i = 0; i < K; ++i) {
            G.template segment<T>(i*T) =
                this->m_function(this->m_time, this->m_state);
        }

        for (size_t iter = 0; iter < m_newton_iters; ++iter) {
            StageVector F = StageVector::Zero();
            StageMatrix J = StageMatrix::Zero();

            for (size_t i = 0; i < K; ++i) {
                VectorT xi = this->m_state;
                for (size_t j = 0; j < K; ++j) {
                    xi += dt * this->m_A(i, j)
                        * G.template segment<T>(j*T);
                }

                double ti = this->m_time + dt * this->m_C(i);
                VectorT fi = this->m_function(ti, xi);

                F.template segment<T>(i*T) =
                    G.template segment<T>(i*T) - fi;

                MatrixT df = numerical_jacobian(ti, xi);

                for (size_t j = 0; j < K; ++j) {
                    MatrixT block = MatrixT::Zero();
                    if (i == j) block.setIdentity();
                    block -= dt * this->m_A(i, j) * df;
                    J.template block<T,T>(i*T, j*T) = block;
                }
            }

            StageVector delta = J.fullPivLu().solve(-F);
            G += delta;

            if (delta.norm() < m_tol) break;
        }

        VectorT high = VectorT::Zero();
        VectorT low  = VectorT::Zero();

        for (size_t i = 0; i < K; ++i) {
            high += this->m_B(i)
                  * G.template segment<T>(i*T);
            low  += m_b_low(i)
                  * G.template segment<T>(i*T);
        }

        VectorT err = dt * (high - low);
        double err_norm = err.norm();

        double factor = m_safety * std::pow(m_tol / (err_norm + 1e-14),
                                            1.0 / m_order_diff);
        factor = std::clamp(factor, 0.2, 5.0);

        if (err_norm <= m_tol) {
            this->m_state += dt * high;
            this->m_time  += dt;
            dt *= factor;
        } else {
            dt *= factor;
        }

        return err_norm;
    }

private:
    Eigen::Matrix<double, K, 1> m_b_low;
    size_t m_newton_iters;
    double m_tol;
    double m_safety;
    static constexpr double m_order_diff = 1.0;

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
