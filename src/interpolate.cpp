#include <interpolate.h>

void interpolate(const Eigen::Vector3d& X, const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d &phi) {
	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3, 2);
	T << X1 - X0, X2 - X0;

	phi.segment<2>(1) = (T.transpose() * T).inverse() * T.transpose() * (X - X0);
	phi(0) = 1 - phi(1) - phi(2);
}