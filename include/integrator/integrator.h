#pragma once
#include <eigen3/Eigen/Dense>
#include <functional>

/**
 * @brief base class for numerical integrators
 * 
 * @tparam T dimension of the function to be integrated
 */
template <class System, size_t T>
class integrator{

protected:
    // maybe dont save this and just pass 
    double m_time;
    Eigen::Matrix<double, T, 1> m_state;
    System m_system;

public:
    integrator(System system): m_system(system){}

    void set_state(const Eigen::Matrix<double,T,1>& state){
        for(int i = 0;i<T;i++)
            m_state(i) = state(i);
    }

    Eigen::Matrix<double, T, 1> get_state(){
        return m_state;
    }


    virtual void step(double dt) = 0;
    //virtual void propagate_for(double time);
};