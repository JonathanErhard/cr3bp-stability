#include "explicit-runge-kutta.h"
#include <fstream>
#include <eigen3/Eigen/Dense>

#define dt 1e-5
#define STEPS 1000000
#define T 2
#define K 3

void rk3_exp(){
    auto exp = [](double t, Eigen::Matrix<double, T, 1> state){
        return state;
    };

    Eigen::Matrix<double, K, K> a = Eigen::Matrix<double, K, K>::Zero();
    a << 0,   0, 0,
         0.5, 0, 0,
         -1,  2, 0;

    Eigen::Matrix<double, K, 1> b;
    b << 1.0/6, 2.0/3, 1.0/6;

    ExplicitRungeKutta<T,K> rk_integrator = ExplicitRungeKutta<T,K>(exp, a, b);

    Eigen::Matrix<double,T,1> state;
    state << 1,2;

    rk_integrator.set_state(state);
    double t = 0.0;

    std::ofstream of;
    of.open("data.csv");

    for(int i=0;i<STEPS;i++){
        t += dt;
        rk_integrator.step(dt);
        state = rk_integrator.get_state();
        of << t;
        for(int i = 0; i< T;i++)
            of << ',' << state[i];
        of << '\n';
    }
}

void euler_exp(){
    Eigen::Matrix<double,T,1> d;
    d << 1.0, 2.0;
    double t = 0;
    for(int i = 0;i<STEPS;i++){
        t += dt;
        d += dt*d;
    }
}

int main(){
    // std::cout << "rk3:\n";
    rk3_exp();
    // std::cout << "euler:\n";
    // euler_exp();
}