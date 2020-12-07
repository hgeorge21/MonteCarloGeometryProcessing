#ifndef GREENS_FUNCTIONS_H
#define GREENS_FUNCTIONS_H

#include <Eigen/Core>

#define PI 3.1415926535897932384626433832795028884L

enum GREEN_FUNC_TYPE {
    HARMONIC,
    YUKAWA
};

// Computes the harmonic or Yukawa potential Green's function
//
// Inputs:
//   X #P by 3 list of query points
//   Y #P by 3 list of points in the ball B_r(X)
//   R #P by 3 list of radius of spheres centered at X
//   c constant used for screened Poisson (not implemented)
//   type Type of Green's function to use
// Output:
//   G Green's function of X, Y are domain of spheres
void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G);

// Computes the harmonic or Yukawa potential Green's function and integral
//
// Inputs:
//   X #P by 3 list of query points
//   Y #P by 3 list of points in the ball B_r(X)
//   R #P by 3 list of radius of spheres centered at X
//   c constant used for screened Poisson (not implemented)
//   type Type of Green's function to use
// Output:
//   G Green's function of X, Y are domain of spheres
//   G_int the integral of Green's function over its domain of spheres
void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G,
        Eigen::VectorXd &G_int);

// Computes the harmonic or Yukawa potential Green's function and integral
//
// Inputs:
//   X #P by 3 list of query points
//   Y #P by 3 list of points in the ball B_r(X)
//   R #P by 3 list of radius of spheres centered at X
//   c constant used for screened Poisson (not implemented)
//   type Type of Green's function to use
// Output:
//   G Green's function of X, Y are domain of spheres
//   G_int the integral of Green's function over its domain of spheres
//   G_grad the gradient of Green's function over its domain of spheres
void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G,
        Eigen::VectorXd &G_int,
        Eigen::MatrixXd &G_grad);

#endif // !GREENS_FUNCTIONS_H
