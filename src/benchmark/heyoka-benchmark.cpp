#include <array>
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>

#include <benchmark/roundtrip_closure.h>
#include <benchmark/surrogate_P1.h>
#include <benchmark/hamiltonian.h>

#include <heyoka/taylor.hpp>
#include <model/cr3bp_heyoka_expressions.h>

#include <dbg.h>

#include <cfg/benchmark-hyper-parameters.h>

using state_type = std::array<double, 6>;
using trajectory_type = std::vector<std::pair<double,state_type>>;

using namespace heyoka;

struct taylor_adaptive_wrappper{
    double m_abs, m_rel;
    taylor_adaptive_wrappper(double abs, double rel){
        m_abs = abs;
        m_rel = rel;
    }

    template<class System>
    trajectory_type integrate(System system, state_type x0, double t0, double tf, double dt){
        std::vector<double> x0_vec;
        for(int i = 0;i<6;i++)
            x0_vec.push_back(x0[i]);
        auto ta = taylor_adaptive<double>(system, x0_vec,kw::tol = INTEGRATOR_PARAMETERS::relative_tolerance); // TODO set tolerance
        ta.set_time(t0);
        auto output = std::get<4>(ta.propagate_for(tf-t0,heyoka::kw::c_output = true));


        trajectory_type traj;

        // if dt is negative, we need to go backwards, which is why tf would be a lower bound instead of an upper bound. I Had a lot of fun figuring this one out ._.
        for (double tm = t0; (dt > 0) ? (tm <= tf) : (tm >= tf); tm += dt) {
            (*output)(tm);
            std::array<double,6> state{
                output->get_output()[0],
                output->get_output()[1], 
                output->get_output()[2],
                output->get_output()[3],
                output->get_output()[4],
                output->get_output()[5]
            };
            traj.emplace_back(tm,state);
        }   

        return traj;
    }
};


int main(){
    // set up integrator + directory to store benchmark results
    const std::string integrator_name = "Heyoka";
    taylor_adaptive_wrappper subject{INTEGRATOR_PARAMETERS::absolute_tolerance,
            INTEGRATOR_PARAMETERS::relative_tolerance};

    const std::string file_path_prefix = "../benchmark_output/" + integrator_name;
    std::filesystem::create_directories(file_path_prefix);
    
    // set up model
    auto cr3bp_model = cr3bp_heyoka_expr(INTEGRATOR_PARAMETERS::mu);
    auto p1_surrogate_model = surrogate_p1_heyoka_expr();
    
    #include<benchmark/run-benchmark.h>
    return 0;
}