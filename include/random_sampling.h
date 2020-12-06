#ifndef RANDOM_SAMPLING_H
#define RANDOM_SAMPLING_H

#include <Eigen/Core>
#include <Eigen/Geometry>

// Randomly samples points from a mesh (V, F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
// Output:
//   P  n by 3 list of query positions

void random_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& P);

#endif // !RANDOM_SAMPLING_H
