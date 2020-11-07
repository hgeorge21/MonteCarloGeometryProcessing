#ifndef WOS_H
#define WOS_H

#include <Eigen/Core>
#include <functional>
#include <igl/AABB.h>

// Computes the Monte Carlo estimate for u(x) using Walk-on-Sphere algorithm
// Solves  Δu = 0 / f at x0.
//
// Inputs:
//   V:     #V x 3 matrix of vertex coordinates
//   F:     #F x 3 matrix of faces, each row is indices into V of a triangle
//   aabb:  AABB tree representation of mesh (V, F)
//   u_hat: estimation function passing in a point, radius of largest ball, and whether on boundary
//   x0:    single point to evaluate u at
// Outputs:
//   u：    value of u at points in x0 using WoS
//
void walk_on_sphere(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const igl::AABB<Eigen::MatrixXd, 3> aabb,
	const std::function<float(Eigen::Vector3d, double, bool)> u_hat, const Eigen::Vector3d& x0, double& u);

#endif // !WOS_H