#include <eigen3/Eigen/Dense>
#include <functional>

/**
 * @brief base class for numerical integrators
 * 
 * @tparam T dimension of the function to be integrated
 */
template <size_t T>
class integrator{

protected:
    // maybe dont save this and just pass 
    double m_time;
    Eigen::Matrix<double, T, 1> m_state;
    std::function<Eigen::Matrix<double, T, 1>(double,Eigen::Matrix<double, T, 1>)> m_function;

public:
    integrator(std::function<Eigen::Matrix<double, T, 1>(double,Eigen::Matrix<double, T, 1>)> function): m_function(function){}

    void set_state(const Eigen::Matrix<double,T,1>& state){
        m_state = state;
    }

    Eigen::Matrix<double, T, 1> get_state(){
        return m_state;
    }


    virtual void step(double dt) = 0;
};