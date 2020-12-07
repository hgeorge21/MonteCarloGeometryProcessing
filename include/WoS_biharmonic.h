#ifndef WOS_BIHARMONIC_H
#define WOS_BIHARMONIC_H

#include <Eigen/Core>

// Solve the following system over space at given points P subject to
// conditions on the given boundary mesh (V,F)
//   ∆^2 u = 0 on Ω
//   u     = g on ∂Ω
//   ∆u    = h on ∂Ω
//
// Inputs:
//   V             #V by 3 list of surface mesh vertex positions
//   F             #F by 3 list of triangles
//   Bg            #V by 1 list of Dirichlet boundary conditions for u
//   Bh            #V by 1 list of Dirichlet boundary conditions for ∆u
//   use_pt_src    whether to use point source for importance sampling
//   point_source  used for importance sampling if enabled
//   P             #P by 3 list of query positions
// Outputs:
//   U             #P by 1 list of values at query positions
void WoS_biharmonic(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& Bg,
        const Eigen::VectorXd& Bh,
        const bool &use_pt_src,
        const Eigen::RowVector3d& point_source,
        const Eigen::MatrixXd& P,
        Eigen::VectorXd& U);

#endif //WOS_BIHARMONIC_H
