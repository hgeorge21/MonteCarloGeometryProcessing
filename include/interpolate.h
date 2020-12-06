#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <Eigen/Core>
#include <Eigen/Geometry>

// Calculate the coefficient of triangle vertices to represent a point in the triangle i.e. linear combination
//
// Inputs:
//   X: point in the triangle
//   X0: first vertex in the triangle
//   X1: second vertex in the triangle
//   X2: third vertex in the triangle
// Output:
//   phi: the coefficients for linear combination

void interpolate(const Eigen::Vector3d& X, const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d& phi);
#endif // !INTERPOLATE_H
