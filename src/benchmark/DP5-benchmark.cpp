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
#include <benchmark/utility.h>
#include <benchmark/surrogate_L1.h>

#include <dbg.h>

#include <cfg/benchmark-hyper-parameters.h>

using state_type = std::array<double, 6>;
namespace odeint = boost::numeric::odeint;
using trajectory_type = std::vector<std::pair<double,state_type>>;

using DP5 = odeint::runge_kutta_dopri5<state_type>;

struct dp5_wrapper{
    odeint::dense_output_runge_kutta<odeint::controlled_runge_kutta<DP5>> m_stepper;
    dp5_wrapper(double abs, double rel){
        m_stepper = make_dense_output(
        abs,
        rel,
        DP5{}
        );
    }

    template<class System>
    trajectory_type integrate(System system, state_type x0, double t0, double tf, double dt){
        trajectory_type traj;
        trajectory_observer<state_type> obs(traj);
        odeint::integrate_const(
            m_stepper,
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
    const std::string integrator_name = "Dormand-Prince-5";
    dp5_wrapper subject{INTEGRATOR_PARAMETERS::absolute_tolerance,
            INTEGRATOR_PARAMETERS::relative_tolerance};

    const std::string file_path_prefix = "../benchmark_output/" + integrator_name;
    std::filesystem::create_directories(file_path_prefix);
    
    // set up model
    std::function<void(state_type&,state_type&,double)> cr3bp_model = [&](const state_type& q,state_type& dq,double t){cr3bp(q,dq,t);};
    std::function<void(const state_type&,state_type&,double)> p1_surrogate_model = [&](const state_type& q,state_type& dq,double t){cr3bp_benchmarks::surrogate_p1_ode(q,dq,t);};

        std::stringstream ss;
    ss << " benchmark results:\n";
    ss << "round-trip closure and hamiltonian:\n";
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
        //auto l2_distances = compare_trajectories_isochronous(fwd_traj,bwd_traj);
        //double last_L2 = l2_distances.back();
        double last_L2 = l2_error(fwd_traj.begin()->second, bwd_traj.back().second);
        
        // calculate hamiltonian error of the forward trajectory
        auto hamiltonian_error = cr3bp_benchmarks::hamiltonian_conservation_benchmark(fwd_traj, INTEGRATOR_PARAMETERS::mu);
        double max_hamiltonian_error = std::max_element(
            hamiltonian_error.begin(),
            hamiltonian_error.end(),
            [](const cr3bp_benchmarks::numeric_error& a, const cr3bp_benchmarks::numeric_error& b){
                return abs(a.rel_error) < abs(b.rel_error);
            }
        )->rel_error;

        ss << INTEGRATOR_PARAMETERS::rbc::names[index] << ": last L2 error = " << last_L2 << "\n";
        ss << INTEGRATOR_PARAMETERS::rbc::names[index] << ": max hamiltonian error (log)= " << -std::log10(max_hamiltonian_error) << "\n";

        // save files for plots
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::rbc::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        save_trajectory(fwd_traj, directory_path + "/fw.csv");
        save_trajectory(bwd_traj, directory_path + "/bw.csv");
        
        save_numeric_error(hamiltonian_error, directory_path + "/hamiltonian_error.csv");
    }

    ss << "\nsurrogate_p1 model benchmark:\n";

    for(int index = 0;index<INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions[index];

        // loop through surrogate_p1 starting positions
        const auto& [exact_traj,estimated_traj] = cr3bp_benchmarks::surrogate_p1_benchmark(
            p1_surrogate_model,
            subject,
            x0,
            0,
            INTEGRATOR_PARAMETERS::integration_time,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
            
        auto l2_distances = compare_trajectories_isochronous(exact_traj,estimated_traj);
        double last_L2 = l2_distances.back().abs_error;

        // safe results to file
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::surrogate_p1::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        ss << INTEGRATOR_PARAMETERS::surrogate_p1::names[index] << ": last error = " << last_L2 << "\n";
        // save trajectories
        save_trajectory(exact_traj, directory_path + "/exact_surrogate.csv");
        save_trajectory(estimated_traj, directory_path + "/estimated_surrogate.csv");

    }
    std::cout << "Files written to " << file_path_prefix << "/<csv-files>.\n";
    std::cout << integrator_name <<  ss.str() << "\n####################\n\n";
    return 0;
}