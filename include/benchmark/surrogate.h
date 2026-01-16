#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <cfg/benchmark-hyper-parameters.h>
#include <functional>

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

namespace cr3bp_benchmarks{

using state_type = std::array<double, 6>;

// vector arithmetics. In hindsight I should have used eigen
std::array<double,3> vec(const state_type& s, int offset) {
    return { s[offset], s[offset+1], s[offset+2] };
}

std::array<double,3> add(const std::array<double,3>& a,
                         const std::array<double,3>& b) {
    return { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
}

std::array<double,3> sub(const std::array<double,3>& a,
                         const std::array<double,3>& b) {
    return { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
}

std::array<double,3> mul(const std::array<double,3>& v, double s) {
    return { v[0]*s, v[1]*s, v[2]*s };
}

double dot(const std::array<double,3>& a,
           const std::array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::array<double,3> cross(const std::array<double,3>& a,
                           const std::array<double,3>& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

double norm(const std::array<double,3>& v) {
    return std::sqrt(dot(v,v));
}

double solveKeplerElliptic(double M, double e) {
    double E = M;
    for (int i = 0; i < 50; ++i)
        E -= (E - e*std::sin(E) - M) / (1.0 - e*std::cos(E));
    return E;
}

double solveKeplerHyperbolic(double M, double e) {
    double H = std::asinh(M / e);
    for (int i = 0; i < 50; ++i)
        H -= (e*std::sinh(H) - H - M) / (e*std::cosh(H) - 1.0);
    return H;
}

// this uses exact keplerian orbits to calculate a reference orbit.
std::vector<std::pair<double, state_type>> surrogate_p1_exact(
    const state_type& initial,
    double mu,
    double duration,
    double grid_resolution
) {
    double samples = duration/grid_resolution;
    auto r0 = vec(initial, 0);
    auto v0 = vec(initial, 3);

    double r0mag = norm(r0);
    double v0mag = norm(v0);

    auto h = cross(r0, v0);
    double hmag = norm(h);

    auto evec = sub(mul(cross(v0, h), 1.0 / mu),
                    mul(r0, 1.0 / r0mag));
    double e = norm(evec);

    double energy = 0.5*v0mag*v0mag - mu / r0mag;
    bool bound = (e < 1.0);

    double a = bound ? -mu / (2.0 * energy)
                     :  mu / (2.0 * std::abs(energy));

    double n = std::sqrt(mu / std::abs(a*a*a));

    auto rhat = mul(evec, 1.0 / e);
    auto qhat = mul(cross(h, rhat), 1.0 / hmag);

    std::vector<std::pair<double, state_type>> trajectory;
    trajectory.reserve(samples + 1);

    for (int i = 0; i <= samples; ++i) {
        double t = duration * i / samples;
        double M = n * t;

        double nu, r;

        if (bound) {
            double E = solveKeplerElliptic(M, e);
            nu = 2.0 * std::atan2(
                std::sqrt(1+e)*std::sin(E/2),
                std::sqrt(1-e)*std::cos(E/2)
            );
            r = a * (1.0 - e*std::cos(E));
        } else {
            double H = solveKeplerHyperbolic(M, e);
            nu = 2.0 * std::atan2(
                std::sqrt(e+1)*std::sinh(H/2),
                std::sqrt(e-1)*std::cosh(H/2)
            );
            r = a * (e*std::cosh(H) - 1.0);
        }

        auto rvec = add(mul(rhat, r*std::cos(nu)),
                        mul(qhat, r*std::sin(nu)));

        double vr = mu / hmag * e * std::sin(nu);
        double vt = mu / hmag * (1.0 + e*std::cos(nu));

        auto vvec = add(mul(rhat, vr),
                        mul(qhat, vt));

        trajectory.push_back({
            t,
            {
                rvec[0], rvec[1], rvec[2],
                vvec[0], vvec[1], vvec[2]
            }
        });
    }

    return trajectory;
}

void surrogate_p1_ode(const state_type& q, state_type& dq, double t){
    const double mu = INTEGRATOR_PARAMETERS::mu;
    // r1 = 0
    const double& x = q[0];
    const double& y = q[1];
    const double& z = q[2];
    const double& dx = q[3];
    const double& dy = q[4];
    const double& dz = q[5];

    double r_abs_1 = 1/sqrt(sqr(x)+sqr(y)+sqr(z)); //absolute value of the position, raised to 2/3

    double ax = - mu * x * r_abs_1* r_abs_1* r_abs_1;
    double ay = - mu * y * r_abs_1* r_abs_1* r_abs_1;
    double az = - mu * z * r_abs_1* r_abs_1* r_abs_1;

    dq[0] = dx;
    dq[1] = dy;
    dq[2] = dz;
    dq[3] = ax;
    dq[4] = ay;
    dq[5] = az;
}

template<class Integrator>
std::pair<std::vector<std::pair<double,state_type>>,std::vector<std::pair<double,state_type>>> surrogate_p1_benchmark(    Integrator integrator_,
    const state_type &x0_s,
    double t0,
    double tf,
    double dt){
    
    using trajectory_type = std::vector<std::pair<double,state_type>>;
    namespace odeint = boost::numeric::odeint;

    state_type x0(x0_s);
    
    auto exact_traj = surrogate_p1_exact(x0, INTEGRATOR_PARAMETERS::mu, tf,dt);
    
    trajectory_type estimated_traj;
    trajectory_observer<state_type> estimation_obs(estimated_traj);
    std::function<void(state_type&,state_type&,double)> system = [&](const state_type& q,state_type& dq,double t){surrogate_p1_ode(q,dq,t);};

    odeint::integrate_const(
        integrator_,
        system,
        x0,
        t0,
        tf,
        dt,
        estimation_obs
    );

    return {exact_traj,estimated_traj};
}

};