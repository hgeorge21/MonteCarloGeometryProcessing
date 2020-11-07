#ifndef FIND_CLOSEST_POINT
#define FIND_CLOSEST_POINT

#include <Eigen/Core>

// Finds the closest point on the boundary of mesh to query point p
//
// Input:
//   V: #V x 3 matrix of vertex coordinates
//   F: #F x 3 matrix of indices into V representing face of triangle
//   p: query point
// Output:
//   d: closest square distance
//   c: closest point's coordinate

void find_closest_point(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::Vector3d& p, double &d, Eigen::Vector3d& c);

#endif // !FIND_CLOSEST_POINT
