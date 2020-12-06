#include <pts_on_spheres.h>
#include <chrono>
#include <random>

void pts_on_spheres(const Eigen::VectorXd& R, Eigen::MatrixXd& X) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist_theta(0., std::_Pi);
	std::uniform_real_distribution<double> dist_phi(0., 2 * std::_Pi);

	auto& sine = [](double x) { return sin(x); };
	auto& cosine = [](double x) { return cos(x); };

	Eigen::VectorXd thetas = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_theta(generator); });
	Eigen::VectorXd phis = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_phi(generator); });

	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(X.rows(), 3);
	res.col(0) = R.cwiseProduct((thetas.unaryExpr(sine)).cwiseProduct(phis.unaryExpr(cosine)));
	res.col(1) = R.cwiseProduct((thetas.unaryExpr(sine)).cwiseProduct(phis.unaryExpr(sine)));
	res.col(2) = R.cwiseProduct((thetas.unaryExpr(cosine)));
	
	X = X + res;
}