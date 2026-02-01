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

    #include<benchmark/run-benchmark.h>
    return 0;
}