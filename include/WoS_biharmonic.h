#ifndef WOS_BIHARMONIC_H
#define WOS_BIHARMONIC_H

#include <Eigen/Core>

// Solve ∆u = 0 over space at given poinst P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles
//   B  #V by 1 list of Dirichlet boundary conditions
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions

void WoS_biharmonic(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& Bh,
        const Eigen::VectorXd& Bg,
        const bool &use_src_pt,
        const Eigen::RowVector3d& source_point,
        const Eigen::MatrixXd& P,
        Eigen::VectorXd& U);

#endif //WOS_BIHARMONIC_H
