#include <iostream>
#include <model/cr3bp_heyoka_expressions.h>
#include <model/CR3BP.h>
#include <eigen3/Eigen/Dense>


using namespace heyoka;

int main()
{
    double mu_val = 1.215e-2;

    auto ode = JE::cr3bp_vel(mu_val);

    // Extract variables and RHS expressions
    std::vector<heyoka::expression> vars;
    std::vector<heyoka::expression> rhs;

    vars.reserve(ode.size());
    rhs.reserve(ode.size());

    for (auto &p : ode) {
        vars.push_back(p.first);
        rhs.push_back(p.second);
    }

    // Build a compiled function f(vars → rhs)
    heyoka::v40::cfunc<double> f(rhs, vars);

    std::vector<double> x0 = {
        8.2339081983651485e-1,// x
        -1.9017764504099543e-28,// y
        9.8941366235910004e-4,// z
        -2.3545391932685812e-15,// vx
        1.2634272983881797e-1,// vy
        2.2367029429442455e-16// vz
    };
    Eigen::Matrix<double,6,1> x0_eigen = Eigen::Map<Eigen::Matrix<double, 6, 1> >(x0.data());
    Eigen::Matrix<double,6,1> derivative;
    CR3BP::calculate_derivative_rotating<double>(x0_eigen,mu_val,derivative);

    std::cout << "state:\n";
    for (auto v : x0)
        std::cout << v << "\n";

    std::array<double,6> out{};

    f(out, x0);

    std::cout << "Derivative my model:\n" << derivative << '\n';

    std::cout << "Derivative heyoka:\n";

    for (auto v : out)
        std::cout << v << "\n";

    return 0;
}