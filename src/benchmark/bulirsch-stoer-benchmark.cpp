#include <array>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>

#include <model/CR3BP.h>
#include <benchmark/roundtrip_closure.h>
#include <benchmark/surrogate_P1.h>
#include <benchmark/hamiltonian.h>

#include <dbg.h>

#include <cfg/benchmark-hyper-parameters.h>

using state_type = std::array<double, 6>;
using trajectory_type = std::vector<std::pair<double,state_type>>;
namespace odeint = boost::numeric::odeint;

using bulirsch_stoer = odeint::bulirsch_stoer_dense_out< state_type >;

struct bulrisch_stoer_wrappper{
    double m_abs, m_rel;
    bulrisch_stoer_wrappper(double abs, double rel){
        m_abs = abs;
        m_rel = rel;
    }

    template<class System>
    trajectory_type integrate(System system, state_type x0, double t0, double tf, double dt){
        bulirsch_stoer stepper = bulirsch_stoer(m_abs,m_rel,0.0 , 0.0);
        trajectory_type traj;
        trajectory_observer<state_type> obs(traj);
        odeint::integrate_const(
            stepper,
            system,
            x0,
            t0,
            tf,
            dt,
            obs
        );
        return obs.traj;
    }
};


int main(){
    // set up integrator + directory to store benchmark results
    const std::string integrator_name = "bulirsch-stoer";
    bulrisch_stoer_wrappper subject{INTEGRATOR_PARAMETERS::absolute_tolerance,
            INTEGRATOR_PARAMETERS::relative_tolerance};

    const std::string file_path_prefix = "../benchmark_output/" + integrator_name;
    std::filesystem::create_directories(file_path_prefix);
    
    // set up models
    std::function<void(const state_type&,state_type&,double)> cr3bp_model = [&](const state_type& q,state_type& dq,double t){cr3bp(q,dq,t);};
    std::function<void(const state_type&,state_type&,double)> p1_surrogate_model = [&](const state_type& q,state_type& dq,double t){cr3bp_benchmarks::surrogate_p1_ode(q,dq,t);};
    
    // loop through round-back-closure starting states
    for(int index = 0;index<INTEGRATOR_PARAMETERS::rbc::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::rbc::starting_positions[index];

        const auto& [fwd_traj,bwd_traj] = cr3bp_benchmarks::roundtrip_closure_benchmark(
            cr3bp_model,
            subject,
            INTEGRATOR_PARAMETERS::rbc::starting_positions[index],
            0.0,
            INTEGRATOR_PARAMETERS::integration_time,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
        double max_l2 = compare_trajectories_isochronic(fwd_traj,bwd_traj);
        
        // calculate hamiltonian error of the forward trajectory
        auto hamiltonian_error = cr3bp_benchmarks::hamiltonian_conservation_benchmark(fwd_traj, INTEGRATOR_PARAMETERS::mu);
        
        // save files for plots
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::rbc::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        save_trajectory(fwd_traj, directory_path + "/fw_1.csv");
        save_trajectory(bwd_traj, directory_path + "/bw_1.csv");
        trajectory_type joined_traj;
        std::copy(fwd_traj.begin(),fwd_traj.end(),std::back_inserter(joined_traj));
        joined_traj.insert(joined_traj.end(),bwd_traj.begin(),bwd_traj.end());
        save_trajectory(joined_traj, directory_path + "/fw_bw.csv");
        
        save_numeric_error(hamiltonian_error, directory_path + "/hamiltonian_error.csv");
    }

    for(int index = 0;index<INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions[index];

        // loop through surrogate_p1 starting positions
        const auto& [exact_traj,estimated_traj] = cr3bp_benchmarks::surrogate_p1_benchmark(
            p1_surrogate_model,
            subject,
            x0,
            0,
            100,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
            
        // safe results to file
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::surrogate_p1::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        // save trajectories
        save_trajectory(exact_traj, directory_path + "/exact_surrogate.csv");
        save_trajectory(estimated_traj, directory_path + "/estimated_surrogate.csv");
        
        // save both trajectories one after the other for plotting purposes
        trajectory_type t1;
        std::copy(exact_traj.begin(),exact_traj.end(),std::back_inserter(t1));
        t1.insert(t1.end(),estimated_traj.begin(),estimated_traj.end());
        save_trajectory(t1, directory_path + "/exact_and_estimated_surrogate.csv");
    }
    std::cout << "Files written to " << file_path_prefix << "/<csv-files>.\n";
    return 0;
}