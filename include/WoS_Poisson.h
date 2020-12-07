#ifndef MONTECARLOPDE_WOS_POISSON_H
#define MONTECARLOPDE_WOS_POISSON_H

#include <Eigen/Core>

void WoS_Poisson(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& B,
        std::function<double(const Eigen::Vector3d &)> f,
        const bool &use_src_pt,
        const Eigen::RowVector3d& source_point,
        const Eigen::MatrixXd& P,
        Eigen::VectorXd& U);

#endif //MONTECARLOPDE_WOS_POISSON_H
