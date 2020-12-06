#include <pts_on_spheres.h>
#include <random>

void pts_on_spheres(const Eigen::VectorXd& R, Eigen::MatrixXd& X) {
	std::default_random_engine generator;
	std::uniform_real_distribution<double> dist_theta(0., PI);
	std::uniform_real_distribution<double> dist_phi(0., 2 * PI);

	const auto& sine = [](double x) { return sin(x); };
	const auto& cosine = [](double x) { return cos(x); };

	Eigen::VectorXd thetas = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_theta(generator); });
	Eigen::VectorXd phis = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_phi(generator); });

	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(X.rows(), 3);
	res.col(0) = R.cwiseProduct((thetas.unaryExpr(sine)).cwiseProduct(phis.unaryExpr(cosine)));
	res.col(1) = R.cwiseProduct((thetas.unaryExpr(sine)).cwiseProduct(phis.unaryExpr(sine)));
	res.col(2) = R.cwiseProduct((thetas.unaryExpr(cosine)));
	
	X = X + res;
}