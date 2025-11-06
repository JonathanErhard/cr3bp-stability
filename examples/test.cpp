#include<sstream>
#include<fstream>
#include<iostream>

#define dbg(x) std::cout << #x << ": " << x << '\n'
#include<model/CR3BP.h>

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
                0.0,
                0.0,
                0.055,
                -0.0;

    Eigen::Matrix<double,6,1> state_L3_backwards;
    state_L3_backwards << 0.923242292247814,0.103708031297204,0.0670944875335092,-0.122382429529801,-0.0484255404818714,-0.031329172693031;

    Eigen::Matrix<double,6,1> state_rnd_backwards;
    state_rnd_backwards << -0.429702389263953,0.678979090260221,0.073379385666045,-0.131232922917919,-0.332455139255467,-0.0357207393086383;

    Eigen::Matrix<double,6,1> state = state_rnd;
    double t = 0.0; 
    
    std::ofstream of;
    of.open("CR3BP.csv");
    of << "t,x,y,z,dx,dy,dz\n"; //csv header for df read in python
    
    for(int i = 0;i<STEPS;i++){
        t+=dt;
        Eigen::Matrix<double,6,1> derivative;
        calculate_derivative_rotating<double>(state,mu,derivative,i%100);
        state += dt*derivative;
        if(i%100 == 0){
            Eigen::Matrix<double,6,1> state_rotating = state;
            of << t << ',' << state_rotating.transpose().format(csv_format) << '\n';
            std::cout << "state: " << state_rotating(0) << "," << state_rotating(1) << "," << state_rotating(3) << "," << state_rotating(4) << '\n';
            std::cout << "deriv: " << derivative(0) << "," << derivative(1) << "," << derivative(3) << "," << derivative(4) << '\n';
            std::cout << calculate_hamiltonian_rotating(state_rotating, mu) << '\n';
        }
    }
    of.close();
}