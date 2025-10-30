#include "integrator.h"

/**
 * @brief base class for (semi) explicit rk and implicit runge kutta methods.
 * 
 * @tparam T dimension of the function, that is being integrated
 * @tparam K dimension of the RK weight matrix (triangular for explicit, any for implicit)
 */

template <size_t T, size_t K>
class RungeKutta : public integrator<T>{
protected:
    const Eigen::Matrix<double, K, K> m_A;
    const Eigen::Matrix<double, K, 1> m_B;
    Eigen::Matrix<double, K, 1> m_C;
public:
    RungeKutta(
        std::function<Eigen::Matrix<double, T, 1>(double t, Eigen::Matrix<double, T, 1>)> function,
        const Eigen::Matrix<double, K, K> a,
        const Eigen::Matrix<double, K, 1> b):
        m_A(a),
        m_B(b),
        integrator<T>(function)
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