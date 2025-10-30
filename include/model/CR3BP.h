#include<eigen3/Eigen/Dense>

namespace CR3BP{
#define sqr(x) ((x)*(x))

/**
 * @brief calculate derivative of the state x,y,z,dx,dy,dz in the rotating frame of a CR3BP
 * 
 * @tparam SCALAR_TYPE 
 * @param state [pos,pos_dot] 
 * @param mu mass parameter [1]
 * @param derivative return [pos_dot,pos_2dot]
 */
template<typename SCALAR_TYPE>
inline void calculate_derivative_rotating(const Eigen::Matrix<SCALAR_TYPE,6,1>& state, const double mu, Eigen::Matrix<SCALAR_TYPE,6,1>& derivative){
    const double mu1 = 1-mu;
    const double mu2 = mu;

    const SCALAR_TYPE& x =  state[0];
    const SCALAR_TYPE& y =  state[1];
    const SCALAR_TYPE& z =  state[2];
    const SCALAR_TYPE& dx = state[3];
    const SCALAR_TYPE& dy = state[4];
    const SCALAR_TYPE& dz = state[5];
    const SCALAR_TYPE r1sqr = sqr(x+mu2)+sqr(y)+sqr(z);
    const SCALAR_TYPE r2sqr = sqr(x-mu1)+sqr(y)+sqr(z);

    const SCALAR_TYPE r1_n3_2 = pow(r1sqr,-1.5);
    const SCALAR_TYPE r2_n3_2 = pow(r2sqr,-1.5);

    // derivative of U_bar wrt the radii
    const double dUr_2 = 1/2*(mu1+mu1*r1_n3_2 + mu2+mu2*r2_n3_2);
    const double dUx = (r1_n3_2-r2_n3_2)*mu1*mu2 + 2*x*dUr_2;
    const double dUy = 2*y*dUr_2;
    const double dUz = 2*z*dUr_2;
    
    derivative[3] =  2*dy - dUx;
    derivative[4] = -2*dx - dUy;
    derivative[5] = -dUz;

    derivative[0] = state[3];
    derivative[1] = state[4];
    derivative[2] = state[5];
}

/**
 * @brief Returns a rotation matrix, that maps from the rotating frame of the CR3BP to the inertial frame
 * 
 * @tparam SCALAR_TYPE 
 * @param t normalized time [T/2PI]
 * @param rotation_matrix return
 */
template<typename SCALAR_TYPE>
inline Eigen::Matrix<SCALAR_TYPE,3,3> r_rotating_inertial(SCALAR_TYPE t){
    Eigen::Matrix<SCALAR_TYPE,3,1> axis(1,0,0);
    return Eigen::AngleAxis<SCALAR_TYPE>(t,axis).toRotationMatrix();
}

/**
 * @brief Returns a rotation matrix, that maps from the inertial frame to the rotating frame of the CR3BP to
 * 
 * @tparam SCALAR_TYPE 
 * @param t normalized time [T/2PI]
 * @param rotation_matrix return
 */
template<typename SCALAR_TYPE>
inline void r_inertial_rotating(SCALAR_TYPE t,Eigen::Matrix<SCALAR_TYPE,3,3>& rotation_matrix){
    Eigen::Matrix<SCALAR_TYPE,3,1> axis(-1,0,0);
    rotation_matrix = Eigen::AngleAxis<SCALAR_TYPE>(t,axis);
}

/**
 * @brief propagates state in a CR3BP system with mass parameter mu for dt time units
 * 
 * @tparam SCALAR_TYPE 
 * @param state [pos,pos_dot]
 * @param mu [1]
 * @param t [T/PI]
 * @param dt [T/PI]
 */
template<typename SCALAR_TYPE>
inline void propagate_CR3BP(Eigen::Matrix<SCALAR_TYPE,6,1>& state_inertial, double mu, double t, SCALAR_TYPE dt){
    Eigen::Matrix<SCALAR_TYPE,6,1> derivative;
    Eigen::Matrix<SCALAR_TYPE,6,6> r_inertial_rotating = Eigen::Matrix<SCALAR_TYPE,6,6>::Zero();
    r_inertial_rotating.block(0,0,3,3) = r_inertial_rotating(t);
    r_inertial_rotating.block(3,3,3,3) = r_inertial_rotating(t);
    auto state_rotating = r_inertial_rotating*state_inertial;
    calculate_derivative_rotating(state_rotating,mu,derivative);
    Eigen::Matrix<SCALAR_TYPE,6,6> r_rotating_inertial = Eigen::Matrix<SCALAR_TYPE,6,6>::Zero();
    r_rotating_inertial.block(0,0,3,3) = r_rotating_inertial(t);
    r_rotating_inertial.block(3,3,3,3) = r_rotating_inertial(t);
    state_rotating = state_rotating + dt*derivative;
    state_inertial = r_rotating_inertial*state_rotating;
}

template<typename SCALAR_TYPE>
inline SCALAR_TYPE calculate_hamiltonian_rotating(Eigen::Matrix<SCALAR_TYPE,6,1>& state_rotating, double mu){
    const auto& state = state_rotating;

    const double mu1 = 1-mu;    
    const double mu2 = mu;

    const SCALAR_TYPE& x =  state[0];
    const SCALAR_TYPE& y =  state[1];
    const SCALAR_TYPE& z =  state[2];
    const SCALAR_TYPE& dx = state[3];
    const SCALAR_TYPE& dy = state[4];
    const SCALAR_TYPE& dz = state[5];

    const SCALAR_TYPE r1 = sqrt(sqr(x+mu2)+sqr(y)+sqr(z));
    const SCALAR_TYPE r2 = sqrt(sqr(x-mu1)+sqr(y)+sqr(z));

    const SCALAR_TYPE U_bar = -0.5*(mu1*sqr(r1)+mu2*sqr(r2)) -mu1/r1 -mu2/r2;
    
    // H calculated via Dynamical systems, p.32
    //const SCALAR_TYPE H = 0.5*(sqr(dx)+sqr(dy)+sqr(dz))+U_bar

    // H calculated via chatgpt
    const SCALAR_TYPE px = dx-y;
    const SCALAR_TYPE py = dy+x;
    const SCALAR_TYPE&pz = dz;

    const SCALAR_TYPE H = 0.5*(sqr(px)+sqr(py)+sqr(pz))+y*px-x*py+U_bar;
    return H;
}

// testing stuff
template<typename SCALAR_TYPE>
inline void calculate_derivative_rotating_c(
    const Eigen::Matrix<SCALAR_TYPE,6,1>& state,
    const double mu,
    Eigen::Matrix<SCALAR_TYPE,6,1>& derivative)
{
    const SCALAR_TYPE& x  = state[0];
    const SCALAR_TYPE& y  = state[1];
    const SCALAR_TYPE& z  = state[2];
    const SCALAR_TYPE& dx = state[3];
    const SCALAR_TYPE& dy = state[4];
    const SCALAR_TYPE& dz = state[5];

    const double mu1 = 1.0 - mu;

    // Distances to the primaries
    const SCALAR_TYPE r1 = sqrt( (x + mu)*(x + mu) + y*y + z*z );
    const SCALAR_TYPE r2 = sqrt( (x - 1 + mu)*(x - 1 + mu) + y*y + z*z );

    const SCALAR_TYPE r1_3 = r1*r1*r1;
    const SCALAR_TYPE r2_3 = r2*r2*r2;

    // Accelerations in rotating frame
    const SCALAR_TYPE ddx = 2.0*dy + x - mu1*(x + mu)/r1_3 - mu*(x - 1 + mu)/r2_3;

    const SCALAR_TYPE ddy = -2.0*dx + y - mu1*y/r1_3 - mu*y/r2_3;

    const SCALAR_TYPE ddz = - mu1*z/r1_3 - mu*z/r2_3;

    // Derivative vector
    derivative[0] = dx;
    derivative[1] = dy;
    derivative[2] = dz;
    derivative[3] = ddx;
    derivative[4] = ddy;
    derivative[5] = ddz;
}

template<typename SCALAR_TYPE>
inline SCALAR_TYPE calculate_hamiltonian_rotating_c(
    const Eigen::Matrix<SCALAR_TYPE,6,1>& state_rotating,
    double mu)
{
    const double mu1 = 1.0 - mu;

    const SCALAR_TYPE& x  = state_rotating[0];
    const SCALAR_TYPE& y  = state_rotating[1];
    const SCALAR_TYPE& z  = state_rotating[2];
    const SCALAR_TYPE& dx = state_rotating[3];
    const SCALAR_TYPE& dy = state_rotating[4];
    const SCALAR_TYPE& dz = state_rotating[5];

    const SCALAR_TYPE r1 = sqrt( (x + mu)*(x + mu) + y*y + z*z );
    const SCALAR_TYPE r2 = sqrt( (x - 1 + mu)*(x - 1 + mu) + y*y + z*z );

    // Canonical momenta
    const SCALAR_TYPE px = dx - y;
    const SCALAR_TYPE py = dy + x;
    const SCALAR_TYPE pz = dz;

    // Hamiltonian
    const SCALAR_TYPE H =
        0.5 * (px*px + py*py + pz*pz)
        + y*px - x*py
        - (mu1 / r1 + mu / r2);

    return H;
}

template<typename SCALAR_TYPE>
inline SCALAR_TYPE calculate_hamiltonian_rotating_c_2(
    const Eigen::Matrix<SCALAR_TYPE,6,1>& state_rotating,
    double mu)
{
    const double mu1 = 1.0 - mu;

    const SCALAR_TYPE& x  = state_rotating[0];
    const SCALAR_TYPE& y  = state_rotating[1];
    const SCALAR_TYPE& z  = state_rotating[2];
    const SCALAR_TYPE& dx = state_rotating[3];
    const SCALAR_TYPE& dy = state_rotating[4];
    const SCALAR_TYPE& dz = state_rotating[5];
    
    const SCALAR_TYPE r1 = sqrt( (x + mu)*(x + mu) + y*y + z*z );
    const SCALAR_TYPE r2 = sqrt( (x - 1 + mu)*(x - 1 + mu) + y*y + z*z );

    const SCALAR_TYPE px = dx - y;
    const SCALAR_TYPE py = dy + x;
    const SCALAR_TYPE pz = dz;

    const SCALAR_TYPE H = 0.5 * (px*px + py*py + pz*pz) + y*px - x*py - (mu1 / r1 + mu / r2);

    return H;
}

#undef sqr
}