#pragma once
#include <array>
#include <vector>
#include<eigen3/Eigen/Dense>

#include <model/CR3BP.h>
#include <benchmark/benchmark_namespace.h>

namespace cr3bp_benchmarks{
    struct numeric_error{
        double time;
        double abs_error;
        double rel_error;
    };
    
    std::vector<numeric_error> hamiltonian_conservation_benchmark(const std::vector<std::pair<double,state_type>>& trajectory, double mu){
        std::vector<numeric_error> result;
        auto& starting_state_arr = trajectory[0].second;
        Eigen::Matrix<double,6,1> starting_state_eigen = (Eigen::VectorXd(6) << 
        starting_state_arr[0],
        starting_state_arr[1],
        starting_state_arr[2],
        starting_state_arr[3],
        starting_state_arr[4],
        starting_state_arr[5]).finished();
        double starting_hamiltonian = CR3BP::calculate_hamiltonian_rotating(starting_state_eigen,mu);

        for(const auto& [time,state_arr] : trajectory){
            Eigen::Matrix<double,6,1> state_eigen = (Eigen::VectorXd(6) << state_arr[0],state_arr[1],state_arr[2],state_arr[3],state_arr[4],state_arr[5]).finished();
            double ham = CR3BP::calculate_hamiltonian_rotating(state_eigen,mu);
            double abs_error = abs(ham - starting_hamiltonian);
            double rel_error = abs(abs_error/starting_hamiltonian);
            result.emplace_back(time,abs_error,rel_error);
        }
        return result;
    }
};