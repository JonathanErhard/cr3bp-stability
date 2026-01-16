// this file was inspired by heyokas benchmark implementation

#include <array>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <model/CR3BP.h>

using state_type = std::array<double, 6>;
namespace odeint = boost::numeric::odeint;
using trajectory_type = std::vector<std::pair<double,state_type>>;

struct trajectory_observer {
    trajectory_type &traj;

    explicit trajectory_observer(trajectory_type &t) : traj(t) {}

    void operator()(const state_type &x, double t) const {
        traj.push_back({t,x});
    }
};


double l2_error(const state_type &a, const state_type &b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(sum);
}

double linf_error(const state_type &a, const state_type &b) {
    double max_err = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    return max_err;
}


void pcr3bp(const state_type &q, state_type &dq, double)
{
    Eigen::Matrix<double,6,1> state;
    for(int i = 0;i<6;i++) state(i) = q[i];
    Eigen::Matrix<double,6,1> derivative;
    CR3BP::calculate_derivative_rotating<double>(state,1.215e-2,derivative);
}

template <typename System>
void forward_backward_error(
    System system,
    const state_type &initial_state,
    double t0,
    double tf,
    double dt,
    double abs_tol,
    double rel_tol
) {
    using error_stepper_type =
        odeint::runge_kutta_fehlberg78<state_type>;

    auto controlled_stepper =
        odeint::make_controlled<error_stepper_type>(abs_tol, rel_tol);

    // ---- Forward integration ----
    state_type x_fwd = initial_state;
    trajectory_type forward_traj;

    odeint::integrate_adaptive(
        controlled_stepper,
        system,
        x_fwd,
        t0,
        tf,
        dt,
        trajectory_observer(forward_traj)
    );

    state_type x_bwd = initial_state;
    x_bwd[2] *= -1.0;
    x_bwd[3] *= -1.0;

    trajectory_type backward_traj;

    odeint::integrate_adaptive(
        controlled_stepper,
        system,
        x_bwd,
        t0,
        tf,
        dt,
        trajectory_observer(backward_traj)
    );

    // ---- Error analysis ----
    std::size_t n = std::min(forward_traj.size(), backward_traj.size());

    double max_linf = 0.0;
    double max_l2 = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        std::cout << "t_forward: " << forward_traj[i].first << ", t_backward: " << backward_traj[i].first << '\n';
        max_linf = std::max(
            max_linf,
            linf_error(forward_traj[i].second, backward_traj[i].second)
        );
        max_l2 = std::max(
            max_l2,
            l2_error(forward_traj[i].second, backward_traj[i].second)
        );
    }

    std::cout << "Simulated for 2000.0 time units.\n"; 
    std::cout << "Max L2 error   : " << max_l2 << "\n";


}


int main(int argc, char *argv[])
{
    state_type ic = { -0.8, 0.0, 0.0, -0.6276410653920694 };
    forward_backward_error(
        pcr3bp,
        ic,
        0.0,
        2000.0,
        1e-8,
        1e-15,
        1e-15
    );
}
