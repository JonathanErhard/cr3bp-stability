#pragma once
#include <benchmark/utility.h>


template <typename System, typename Integrator, typename state_type = std::array<double,6>>
[[nodiscard]] std::pair<std::vector<std::pair<double,state_type>>,std::vector<std::pair<double,state_type>>> roundtrip_closure_benchmark(
    System system,
    Integrator integrator_,
    const state_type &x01,
    double t0,
    double tf,
    double dt
) {
    using trajectory_type = std::vector<std::pair<double,state_type>>;
    namespace odeint = boost::numeric::odeint;
    // integrate forward
    state_type x0(x01); //why did I do that? I cant remember but im afraid to change it now
    trajectory_type fwd_traj;
    trajectory_observer<state_type> fwd_obs(fwd_traj);

    odeint::integrate_const(
      integrator_,
      system,
      x0,
      t0,
      tf,
      dt,
      fwd_obs
    );
    //trajectory is initialized with the starting state -> cannot be empty so back() is sure to be defined behaviour
    state_type x0_bwd = (fwd_obs.traj.back()).second;

    // integrate backwards
    trajectory_type bwd_traj;
    trajectory_observer<state_type> bwd_obs(bwd_traj);

    // TODO THIS IS WRONG!!!!
    odeint::integrate_const(
      integrator_,
      system,
      x0_bwd,
      tf,
      t0,
      -dt,
      bwd_obs
    );
    return {fwd_obs.traj,bwd_obs.traj};
}