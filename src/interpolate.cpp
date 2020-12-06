#include <interpolate.h>

void interpolate(const Eigen::Vector3d& X, const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d &phi) {
	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3, 2);
	T << X1 - X0, X2 - X0;

	Eigen::Matrix2d tmp = T.transpose() * T;
	phi.segment<2>(1) = tmp.inverse() * T.transpose() * (X - X0);
	phi(0) = 1 - phi(1) - phi(2);
}