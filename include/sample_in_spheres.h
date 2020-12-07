#ifndef SAMPLE_IN_SPHERES_H
#define SAMPLE_IN_SPHERES_H

#define PI 3.1415926535897932384626433832795028884L

#include <Eigen/Core>

// Samples points in a sphere centered at X
//
// Inputs:
//   X  #P by 3 list of points as center of spheres
//   R  #P by 1 list of sphere radius to walk
// Outputs:
//   Y  #P by 3 list of sampled points in the spheres centered at X
void sample_in_spheres(
        const Eigen::MatrixXd& X,
        const Eigen::VectorXd& R,
        Eigen::MatrixXd& Y);

// Samples points in a sphere centered at X with importance sampling based on source point
//
// Inputs:
//   X  #P by 3 list of points as center of spheres
//   R  #P by 1 list of sphere radius to walk
//   point_source point source used for importance sampling
// Outputs:
//   Y  #P by 3 list of sampled points in the spheres centered at X with importance sampling
void sample_in_spheres(
        const Eigen::MatrixXd& X,
        const Eigen::VectorXd& R,
        const Eigen::RowVector3d &point_source,
        Eigen::MatrixXd& Y);

#endif //SAMPLE_IN_SPHERES_H
