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

    #include<benchmark/run-benchmark.h>
    return 0;
}