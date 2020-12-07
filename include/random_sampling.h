#ifndef RANDOM_SAMPLING_H
#define RANDOM_SAMPLING_H

#include <Eigen/Core>
#include <Eigen/Geometry>

// Randomly samples points from a mesh (V, F)
// Adapted from Libigl Fast Winding Number tutorial
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   n  number of initial samples to generate
// Output:
//   P  n by 3 list of query positions
void random_sampling(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const int &n,
        Eigen::MatrixXd& P);

#endif // !RANDOM_SAMPLING_H
