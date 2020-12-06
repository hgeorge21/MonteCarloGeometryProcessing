#include <pts_on_spheres.h>
#include <random>

void pts_on_spheres(const Eigen::VectorXd& R, Eigen::MatrixXd& X) {
	std::default_random_engine generator;
	std::uniform_real_distribution<double> dist_theta(0., std::_Pi);
	std::uniform_real_distribution<double> dist_phi(0., 2 * std::_Pi);

	for (int k = 0; k < X.rows(); k++) {
		Eigen::RowVector3d vec;
		double theta = dist_theta(generator);
		double phi = dist_phi(generator);

		vec << sin(theta) * cos(phi), sin(theta)* sin(phi), cos(theta);
		vec *= R(k);

		X.row(k) = X.row(k) + vec;
	}
}