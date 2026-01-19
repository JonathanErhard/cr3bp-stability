#pragma once
#include <benchmark/utility.h>

namespace cr3bp_benchmarks{
template <typename System, typename Integrator>
std::pair<std::vector<std::pair<double,state_type>>,std::vector<std::pair<double,state_type>>> roundtrip_closure_benchmark(
    System system,
    Integrator integrator_,
    const state_type &x0,
    double t0,
    double tf,
    double dt
) {
    // integrate forward
    state_type x0_fwd(x0); // copy x0 because of const (and rename)
    trajectory_type fwd_traj = integrator_.integrate(
      system,
      x0_fwd,
      t0,
      tf,
      dt
    );

    // integrate backwards
    state_type x0_bwd = (fwd_traj.back()).second;
    trajectory_type bwd_traj = integrator_.integrate(
      system,
      x0_bwd,
      tf,
      t0,
      -dt
    );
    return {fwd_traj,bwd_traj};
}
}; // namespace cr3bp_benchmarks