#pragma once

#include <tuple>
#include <utility>
#include <vector>
#include <cmath>
#include <variant>

#include <heyoka/expression.hpp>
#include <heyoka/config.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/number.hpp>

#include <cfg/benchmark-hyper-parameters.h>

std::vector<std::pair<heyoka::expression,heyoka::expression>> cr3bp_heyoka_expr(double mu_val){
    auto [x, y, z, dx, dy, dz] = heyoka::make_vars("x", "y", "z", "dx", "dy", "dz");
    
    const auto mu1 = 1.-mu_val;
    const auto mu2 = mu_val;
    
    const auto r1sqr = (x+mu2)*(x+mu2) + y*y + z*z;
    const auto r2sqr = (x-mu1)*(x-mu1) + y*y + z*z;
    
    //const auto r1_n3_2 = pow(r1sqr,3./2);
    //const auto r2_n3_2 = pow(r2sqr,3./2);
    const auto r1_n3_2 = r1sqr * sqrt(r1sqr);
    const auto r2_n3_2 = r2sqr * sqrt(r2sqr);
            
    auto xdot = dx;
    auto ydot = dy;
    auto zdot = dz;
    auto dxdot = x - (mu1 * (x + mu2) / r1_n3_2) - (mu2 * (x - mu1) / r2_n3_2) + 2. * dy;
    auto dydot = y - (mu1 * y / r1_n3_2) - (mu2 * y / r2_n3_2) - 2. * dx;
    auto dzdot = - (mu1 * z /r1_n3_2) - (mu2 * z / r2_n3_2);
    
    return {
        prime(x) = xdot,
        prime(y) = ydot,
        prime(z) = zdot,
        prime(dx)= dxdot,
        prime(dy)= dydot,
        prime(dz)= dzdot
    };
}
std::vector<std::pair<heyoka::expression,heyoka::expression>> surrogate_p1_heyoka_expr(){
    const double mu = INTEGRATOR_PARAMETERS::mu;

    auto [x, y, z, dx, dy, dz] = heyoka::make_vars("x", "y", "z", "dx", "dy", "dz");

    const auto r_abs_1 = 1./sqrt(x*x + y*y + z*z); //absolute value of the position, raised to 2/3

    const auto ax = - mu * x * r_abs_1* r_abs_1* r_abs_1;
    const auto ay = - mu * y * r_abs_1* r_abs_1* r_abs_1;
    const auto az = - mu * z * r_abs_1* r_abs_1* r_abs_1;

    return{
        prime(x) = dx,
        prime(y) = dy,
        prime(z) = dz,
        prime(dx)= ax,
        prime(dy)= ay,
        prime(dz)= az
    };
}