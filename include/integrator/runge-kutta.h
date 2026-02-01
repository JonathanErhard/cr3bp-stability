#pragma once
#include "integrator.h"

/**
 * @brief base class for (semi) explicit rk and implicit runge kutta methods.
 * 
 * @tparam T dimension of the function, that is being integrated
 * @tparam K dimension of the RK weight matrix (triangular for explicit, any for implicit)
 */

template <class System,size_t T, size_t K>
class RungeKutta : public integrator<System,T>{
protected:
    const Eigen::Matrix<double, K, K> m_A;
    const Eigen::Matrix<double, K, 1> m_B;
    Eigen::Matrix<double, K, 1> m_C;
public:
    RungeKutta(
        System system,
        const Eigen::Matrix<double, K, K> a,
        const Eigen::Matrix<double, K, 1> b):
        m_A(a),
        m_B(b),
        integrator<System,T>(system)
        {
            // set c weights
            for(int i = 0;i<K;i++){
                m_C(i) = 0;
                for(int j = 0;j<K;j++){
                    m_C(i) += m_A(i,j);
                }
            }
        }
};