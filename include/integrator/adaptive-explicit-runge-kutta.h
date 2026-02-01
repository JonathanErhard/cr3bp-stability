#pragma once
#include <array>
#include <Eigen/Dense>

#include "runge-kutta.h"

template <class System, size_t T, size_t K>
class AdaptiveExplicitRungeKutta : public RungeKutta<System, T, K> {
public:
    using VectorT = Eigen::Matrix<double, T, 1>;

    AdaptiveExplicitRungeKutta(
        System system,
        const Eigen::Matrix<double, K, K>& a,
        const Eigen::Matrix<double, K, 1>& b_high,
        const Eigen::Matrix<double, K, 1>& b_low,
        double tol = 1e-8,
        double safety = 0.9
    )
        : RungeKutta<System, T, K>(system, a, b_high),
          m_b_low(b_low),
          m_tol(tol),
          m_safety(safety)
    {
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = i; j < K; ++j) {
                if (a(i, j) != 0.0) {
                    std::cerr << "Explicit RK requires lower-triangular A\n";
                    std::exit(-1);
                }
            }
        }
    }

    double step(double& dt) {
        std::array<VectorT, K> g;

        for (size_t i = 0; i < K; ++i) {
            VectorT input = this->m_state;
            for (size_t j = 0; j < i; ++j) {
                input += dt * this->m_A(i, j) * g[j];
            }
            g[i] = this->m_system(this->m_time + dt * this->m_C(i), input);
        }

        VectorT high = VectorT::Zero();
        VectorT low  = VectorT::Zero();

        for (size_t i = 0; i < K; ++i) {
            high += this->m_B(i) * g[i];
            low  += m_b_low(i) * g[i];
        }

        VectorT err = dt * (high - low);
        double err_norm = err.norm();

        // timestep adaptation
        double factor = m_safety * std::pow(m_tol / (err_norm + 1e-14),
                                            1.0 / m_order_diff);
        factor = std::clamp(factor, 0.2, 5.0);

        if (err_norm <= m_tol) {
            this->m_state += dt * high;
            this->m_time  += dt;
            dt *= factor;
            return err_norm;
        } else {
            dt *= factor;
            return err_norm;
        }
    }

private:
    Eigen::Matrix<double, K, 1> m_b_low;
    double m_tol;
    double m_safety;
    static constexpr double m_order_diff = 1.0;
};
