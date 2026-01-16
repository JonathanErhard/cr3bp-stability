#include <iostream>
#include <vector>

#include <heyoka/taylor.hpp>
#include <model/cr3bp_heyoka_expressions.h>

using namespace heyoka;

int main()
{

    const double GRIDSIZE = 0.01;

    double mu_val = 1.215e-2;

    auto eqs = JE::cr3bp_vel(mu_val);
		
    //id 940 halo
//    std::vector<double> x0 = {
//        8.9328742560088747E-1,
//        4.3207679782643764E-27,
//        1.9803451083236487E-1,
//        -2.8045785366944001E-14,
//        1.9828989893836621E-1,
//        1.0076432096249845E-13
//    };

    //id 1147 halo
    std::vector<double> x0 = {
        0,-mu_val, 0, 
        0, 1e-1,1e-1
    };

    auto ta = taylor_adaptive<double>{ eqs, x0,kw::tol = 1e-8};

    double tf = 0.1;//2.7430007981241529E+0;

    auto output = std::get<4>(ta.propagate_for(tf,kw::c_output = true));

    std::cout << "t,x,y,z,dx,dy,dz\n"; 
    for (double tm = 0; tm<tf; tm+=GRIDSIZE) {

        (*output)(tm);
        // Print it out:
        std::cout << tm << ", "
        << output->get_output()[0] << ", "
        << output->get_output()[1] << ", " 
        << output->get_output()[2] << ", "
        << output->get_output()[3] << ", "
        << output->get_output()[4] << ", "
        << output->get_output()[5] << '\n';
    }

    return 0;
}