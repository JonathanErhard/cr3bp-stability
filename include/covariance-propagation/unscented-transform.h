#ifndef _UT_
#define _UT_

/**
 * The ukf equations used for the implementation are described in [1] https://groups.seas.harvard.edu/courses/cs281/papers/unscented.pdf
 * for more information discribing the code I refer to https://news.ycombinator.com/item?id=14339899
 */
#include <eigen3/Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <array>
#include <iostream>
#include <string>


#define STATE_DIM 13
#define MEASUREMENT_DIM 7

// UKF constants
const double alpha = 1e-3;  // spead of sigma points
const double kappa = 0;     // TODO
const double beta = 2;      // parameter discribing the distribution. For gaus 2 is optimal usually used [1]

const double lambda = alpha * alpha * (STATE_DIM + kappa) - STATE_DIM;
const double gamma = (STATE_DIM + lambda) / 2;

std::array<Eigen::Matrix<SCALAR_TYPE,STATE_DIM,1>,2*STATE_DIM+1> sigma_points;

Eigen::Matrix<SCALAR_TYPE, STATE_DIM, STATE_DIM> P;                              // state covariance matrix
Eigen::Matrix<SCALAR_TYPE, STATE_DIM, STATE_DIM> Q_process_noise_cov;            // process nose
Eigen::Matrix<SCALAR_TYPE, MEASUREMENT_DIM, MEASUREMENT_DIM> R_meas_noise_cov;   // measurement noise

std::array<SCALAR_TYPE,2*STATE_DIM+1> W_m; // weights for mean
std::array<SCALAR_TYPE,2*STATE_DIM+1> W_c; // weights for covariances

// diagonal of the astrobee's moment of inertia
const Eigen::Vector3d I = Eigen::Vector3d(0.153, 0.143, 0.162);

// methods for computing the change in omega using runge kutta 4
Eigen::Matrix<SCALAR_TYPE,3,1> rk4OmegaUpdate(const Eigen::Matrix<SCALAR_TYPE,3,1>& omega, double dt);
Eigen::Matrix<SCALAR_TYPE,3,1> computeOmegaDot(const Eigen::Matrix<SCALAR_TYPE,3,1>& omega);

// recomputes sigmapoints for the current state covariance. They do not get propagated yet!
void inline recompute_sigma_points();

// estimated state vector
Eigen::Matrix<SCALAR_TYPE,STATE_DIM,1> x_hat;  
// measurement vector
Eigen::Matrix<SCALAR_TYPE,MEASUREMENT_DIM,1> z_meas;
// inertia parameters


// calculates the value of point at Now() + dt (dynamics developed by Abhi)
inline void propagate_point(Eigen::Matrix<SCALAR_TYPE, STATE_DIM, 1>& point, double dt);

// predicts the state at Now() + dt and updates the state covariance P accordingly
void predict(double dt);
#endif
