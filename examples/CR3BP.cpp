#include<sstream>
#include<fstream>
#include<iostream>

#define dbg(x) std::cout << #x << ": " << x << '\n'
#include<model/CR3BP.h>
#include<integrator/explicit-runge-kutta.h>

#define dt 1e-5
#define STEPS 1000000
#define T 6
#define K 3

using namespace CR3BP;

int main(){
    const Eigen::IOFormat csv_format(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");

    // setting up RK
    Eigen::Matrix<double, K, K> a;
    a << 0,   0, 0,
            0.5, 0, 0,
            -1,  2, 0;
    Eigen::Matrix<double, K, 1> b;
    b << 1.0/6, 2.0/3, 1.0/6;

    const double mu = 1.215e-2;
    const double L = 3.850e5;
    const double V = 1.025;
    const double PERIOD = 2.361e6;

    std::function<Eigen::Matrix<double, T, 1>(double,Eigen::Matrix<double, T, 1>)> propagate = [=] (double t, Eigen::Matrix<double,T,1> state){
        Eigen::Matrix<double, T, 1> derivative;
        calculate_derivative_rotating_c(state,mu,derivative);
        return derivative;
    };

    ExplicitRungeKutta<T,K> rk_integrator = ExplicitRungeKutta<T,K>(propagate,a,b);

    
    // initial state
    Eigen::Matrix<double,6,1> state_L3;
    state_L3 <<    -0.8269,
                0.0,
                0.0,
                0.0,
                0.055,
                0.0;
    Eigen::Matrix<double,6,1> state_rnd;
    state_rnd <<    -0.8269,
                0.3,
                0.2,
                0.0,
                0.055,
                -0.3;

    Eigen::Matrix<double,6,1> state_L3_backwards;
    state_L3_backwards << 0.923242292247814,0.103708031297204,0.0670944875335092,-0.122382429529801,-0.0484255404818714,-0.031329172693031;

    Eigen::Matrix<double,6,1> state_rnd_backwards;
    state_rnd_backwards << -0.429702389263953,0.678979090260221,0.073379385666045,-0.131232922917919,-0.332455139255467,-0.0357207393086383;

    Eigen::Matrix<double,6,1> state = state_rnd_backwards;
    double t = 0.0; 
    rk_integrator.set_state(state);
    
    std::ofstream of;
    of.open("CR3BP.csv");
    of << "t,x,y,z,dx,dy,dz\n"; //csv header for df read in python
    
    for(int i = 0;i<STEPS;i++){
        t+=dt;
        rk_integrator.step(dt);
        if(i%100 == 0){
            Eigen::Matrix<double,6,1> state_rotating = rk_integrator.get_state();
            Eigen::Matrix<double,6,6> R_rotating_inertial = Eigen::Matrix<double,6,6>::Zero();
            R_rotating_inertial.block(0,0,3,3) = r_rotating_inertial(t);
            R_rotating_inertial.block(3,3,3,3) = r_rotating_inertial(t);
            Eigen::Matrix<double,6,1> state = R_rotating_inertial*state_rotating;
            of << t << ',' << state.transpose().format(csv_format) << '\n';
        }
    }
    of.close();
}