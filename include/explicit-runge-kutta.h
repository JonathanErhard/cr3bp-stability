#include <iostream>

#include "runge-kutta.h"

/**
 * @brief base class for (semi) explicit rk and implicit runge kutta methods.
 * 
 * @tparam T dimension of the function, that is being integrated
 * @tparam K dimension of the RK weight matrix (triangular for explicit, any for implicit)
 */

template <size_t T,size_t K>
class ExplicitRungeKutta : public RungeKutta<T, K>{
public:
    ExplicitRungeKutta(
        std::function<Eigen::Matrix<double, T, 1>(double t, Eigen::Matrix<double, T, 1>)> function,
        const Eigen::Matrix<double, K, K> a,
        const Eigen::Matrix<double, K, 1> b):
        RungeKutta<T,K>(function,a,b)
        {
            for(int i = 0; i<T;i++){
                for(int j = i; j<T;j++){
                    if(a(i,j)!=0.0){    
                        std::cout << "\"a\" must be a lower triangular matrix with diagonal equal to 0!\n";
                        exit(-1);
                    }
                }
            }
        }
    void step(double dt){
        // sometimes this is also called K
        std::array<Eigen::Matrix<double,T,1>,K> g;
        for(int i = 0;i<K;i++){
            Eigen::Matrix<double,T,1> input_i = this->m_state;
            for(int j = 1; j<i;j++){
                input_i += dt*this->m_A(i,j)*g[j];
            }
            g[i] = this->m_function(this->m_time + dt*this->m_C(i),input_i);
        }
        Eigen::Matrix<double, T, 1> derivative = Eigen::Matrix<double, T, 1>::Zero();
        for(int i = 0; i<K;i++){
            derivative += this->m_B(i)*g[i];
        }
        this->m_state += dt*derivative;
        this->m_time  += dt;
    }
};