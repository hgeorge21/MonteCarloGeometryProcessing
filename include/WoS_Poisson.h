#ifndef WOS_POISSON_H
#define WOS_POISSON_H

#include <Eigen/Core>

// Solve the following system over space at given points P subject to
// conditions on the given boundary mesh (V,F)
//   ∆u = f on Ω
//   u  = g on ∂Ω
//
// Inputs:
//   V             #V by 3 list of surface mesh vertex positions
//   F             #F by 3 list of triangles
//   B             #V by 1 list of Dirichlet boundary conditions
//   f             the source term function
//   use_pt_src    whether to use point source for importance sampling
//   point_source  used for importance sampling if enabled
//   P             #P by 3 list of query positions
// Outputs:
//   U             #P by 1 list of values at query positions
void WoS_Poisson(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& B,
        std::function<double(const Eigen::Vector3d &)> f,
        const bool &use_pt_src,
        const Eigen::RowVector3d& point_source,
        const Eigen::MatrixXd& P,
        Eigen::VectorXd& U);

#endif //WOS_POISSON_H
