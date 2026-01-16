#pragma once
#include <vector>
#include <array>
#include <string>
namespace INTEGRATOR_PARAMETERS{
    const double mu = 1.215e-2;
    const double integration_time = 5.0;
    const double grid_resolution = 0.01;
    const double absolute_tolerance = 1.0e-15;
    const double relative_tolerance = 1.0e-15;
    
    namespace rbc{
        const std::vector<std::string> names {"halo_stable","halo_unstable","gateway_L2_southern","gateway_L2_northern"};
        const std::vector<std::array<double, 6>> starting_positions{
            std::array<double,6>{8.2339081983651485E-1,-1.9017764504099543E-28,9.8941366235910004E-4,-2.3545391932685812E-15,1.2634272983881797E-1,2.2367029429442455E-16},
            std::array<double,6>{-4.1456184803140111E-1,2.7726895068890510E-23,9.0753120433295065E-1,-1.1555581216266303E-12,1.4076145460136695E+0,3.9979936124406781E-13},
            std::array<double,6>{1.013395,0,-0.175,0,-0.083,0},
            std::array<double,6>{1.0133950,0.0000000,0.1750000,0.0000000,0.0830000,0.0000000}
    };
    }
    namespace surrogate_p1{
        const std::vector<std::string> names {"orbit_right","orbit_above"};
        std::vector<std::array<double, 6>> starting_positions{
            std::array<double,6>{0.2,0,0,0,0.0001,0},
            std::array<double,6>{0,0.2,0,0.0001,0,0}
        };
    };
}
