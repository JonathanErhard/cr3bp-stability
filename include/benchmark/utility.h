#pragma once
#include <array>
#include <vector>

#include <dbg.h>
#include <model/CR3BP.h>
#include <cfg/benchmark-hyper-parameters.h>
#include <benchmark/hamiltonian.h>

/**
 * @brief observer to store a trajectory created from the boost library
 * 
 * @tparam state_type std::array of some length
 */
template <class state_type = std::array<double,6>>
struct trajectory_observer {
    using trajectory_type = std::vector<std::pair<double,state_type>>;
    trajectory_type &traj;
    
    explicit trajectory_observer(trajectory_type &t) : traj(t) {}
    
    // this operator get called by the integrator to output a state
    void operator()(const state_type &x, double t) const {
        traj.push_back({t,x});
    }
};

/**
 * @brief calculates the L2 error of two states
 * 
 * @tparam state_type
 * @param a state a
 * @param b state b
 * @return double euclidian distance/L2 error of the two states
 */
template <class state_type = std::array<double,6>>
double l2_error(const state_type &a, const state_type &b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(sum);
}

/**
 * @brief calculates the L_inf error of two states
 * 
 * @tparam state_type
 * @param a state a
 * @param b state b
 * @return double l_inf error of the two states
 */
template <class state_type = std::array<double,6>>
double linf_error(const state_type &a, const state_type &b) {
    double max_err = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    return max_err;
}

/**
 * @brief wrapper for the model, so that it is compatible with the boost::odeint steppers. it is somewhat inefficient, but executiontime is not what I am testing.
 * 
 * @param q state
 * @param dq derivative
 */
template<typename state_type = std::array<double,6>>
void cr3bp(const state_type &q, state_type &dq, double)
{
    Eigen::Matrix<double,6,1> state;
    for(int i = 0;i<6;i++) state(i) = q[i];
    Eigen::Matrix<double,6,1> derivative;
    CR3BP::calculate_derivative_rotating<double>(state,INTEGRATOR_PARAMETERS::mu,derivative);
    for(int i = 0;i<6;i++) dq[i] = derivative(i);
}

#include<string>
#include<sstream>

template<typename state_type = std::array<double,6>>
inline std::string trajectory_to_str(const std::vector<std::pair<double,state_type>>& traj){
    std::stringstream ss;
    ss << "t,x,y,z,dx,dy,dz\n";
    for(auto& [time, state]:traj){ 
        ss << time << ", "
        << state[0] << ", "
        << state[1] << ", " 
        << state[2] << ", "
        << state[3] << ", "
        << state[4] << ", "
        << state[5] << '\n';
    }
    return ss.str();
}

#include <fstream>

template<typename state_type = std::array<double,6>>
inline void save_trajectory(const std::vector<std::pair<double,state_type>>& traj,const std::string filename){
    std::ofstream of(filename);
    of<<trajectory_to_str(traj);
    of.close();
}

inline std::string numeric_error_to_str(const std::vector<cr3bp_benchmarks::numeric_error>& errors){
    std::stringstream ss;
    ss << "t,abs,rel\n";
    for(auto& [time, rel, abs]:errors) ss << time << ", " << rel << ", " << abs << '\n'; // could use fmt maybe?
    return ss.str();
}

inline void save_numeric_error(const std::vector<cr3bp_benchmarks::numeric_error>& errors,const std::string filename){
    std::ofstream of(filename);
    of<<numeric_error_to_str(errors);
    of.close();
}

template<class state_type = std::array<double,6>>
double compare_trajectories(
    const std::vector<std::pair<double,state_type>>& reference_trajectory,
    const std::vector<std::pair<double,state_type>>& predicted_trajectory){
    std::size_t n = std::min(reference_trajectory.size(), predicted_trajectory.size()); 

    double max_l2 = 0.0;

    for (std::size_t i = 0; i < n; ++i)
        max_l2 = std::max(
            max_l2,
            l2_error(reference_trajectory[i].second, predicted_trajectory[n-i-1].second)
        );

    return max_l2;
}
