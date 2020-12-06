#ifndef POS_H
#define POS_H

#define PI 3.1415926535897932384626433832795028884L

#include <Eigen/Core>

// Calculates the final position on spheres given radius
//
// Inputs:
//   R  #P by 1 list of sphere radius to walk
// Outputs:
//   X  #P by 3 list of points on sphere to go to 

void sample_on_spheres(
	const Eigen::VectorXd& R,
	Eigen::MatrixXd& X);

#endif // !POS_H