#pragma once
#include <vector>
#include <array>
#include <string>
namespace INTEGRATOR_PARAMETERS{
    const double mu = 1.215e-2;
    const double integration_time = 5.0;
    const double grid_resolution = 0.01;
    const double absolute_tolerance = 1.0e-6;
    const double relative_tolerance = 1.0e-6;
    
    namespace rbc{
        const std::vector<std::string> names {"halo_unstable","halo_stable","gateway_L2_southern","gateway_L2_northern","distant_retrograde","transfer_orbit","near_p1"};
        const std::vector<std::array<double, 6>> starting_positions{
            // normal benchmark
            std::array<double,6>{8.2339081983651485E-1,-1.9017764504099543E-28,9.8941366235910004E-4,-2.3545391932685812E-15,1.2634272983881797E-1,2.2367029429442455E-16},
            std::array<double,6>{-4.1456184803140111E-1,2.7726895068890510E-23,9.0753120433295065E-1,-1.1555581216266303E-12,1.4076145460136695E+0,3.9979936124406781E-13},
            std::array<double,6>{1.0133950,0.0000000,0.1750000,0.0000000,0.0830000,0.0000000},
            std::array<double,6>{1.013395,0,-0.175,0,-0.083,0},
            std::array<double,6>{9.8053865097318893E-1,3.6951257762665743E-25,1.5385698388620773E-22,7.3118271005529372E-11,1.2965664912647605E+0,-3.0869020607535532E-19},
            std::array<double,6>{0.03785,0,0,0,-6.35,0}
            //std::array<double,6>{-0.00183592, -0.012712, -7.00826e-05, 1.56495, -1.8083, 0.0947313}
            // fehlberg test
            /*std::array<double,6>{8.2361338183346766E-1,1.2451399523691531E-28,4.2209389835228948E-2,2.0459797529950219E-15,1.5146779732220547E-1,9.6566104380619562E-15},
            std::array<double,6>{8.2338983102485563E-1,8.2680125862110327E-27,3.6867260822713629E-3,-1.8663094559298700E-15,1.2655277462611700E-1,9.0845361544864190E-16},
            std::array<double,6>{8.2338895459048989E-1,-8.1744493465225147E-27,5.0329913747852754E-3,-1.1681655681695368E-15,1.2674774099650740E-1,1.5972995951622074E-15},
            std::array<double,6>{8.2338787303774807E-1,5.9521807372128352E-29,6.3768259805378012E-3,-1.5930354102885625E-15,1.2700163896647410E-1,8.8947025903502806E-16},
            std::array<double,6>{8.2338663761433939E-1,-5.7850809342153700E-27,7.7176461245161143E-3,-1.0326666270228556E-15,1.2731340401147301E-1,1.8079701862210234E-15}*/
    };
    }
    namespace surrogate_p1{
        const std::vector<std::string> names {"orbit_right","orbit_above"};
        std::vector<std::array<double, 6>> starting_positions{
            std::array<double,6>{0.2,0  ,0,0   ,0.01,0},
            std::array<double,6>{0  ,0.2,0,0.001,0   ,0}
        };
    };
}
