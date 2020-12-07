#ifndef SAMPLE_ON_SPHERE_H
#define SAMPLE_ON_SPHERE_H

#include <Eigen/Core>

#define PI 3.1415926535897932384626433832795028884L

// Calculates the final position on spheres given radius
//
// Inputs:
//   R  #P by 1 list of sphere radius
//   X  #P by 3 list of initial center of spheres
// Outputs:
//   X  #P by 3 list of resulting center of spheres
void sample_on_spheres(
	const Eigen::VectorXd& R,
	Eigen::MatrixXd& X);

#endif // !SAMPLE_ON_SPHERE_H