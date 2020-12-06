#include <pts_on_spheres.h>
#include <chrono>
#include <random>

void pts_on_spheres(const Eigen::VectorXd& R, Eigen::MatrixXd& X) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

	std::uniform_real_distribution<double> dist_theta(0., 2 * PI);
	std::uniform_real_distribution<double> dist_u(-1., 1.);

	const auto& sine = [](double x) { return sin(x); };
	const auto& cosine = [](double x) { return cos(x); };
	const auto& smthng = [](double x) { return sqrt(1-x*x); };

	Eigen::VectorXd theta = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_theta(generator); });
	Eigen::VectorXd u = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return dist_u(generator); });

	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(X.rows(), 3);
	res.col(0) = R.cwiseProduct((u.unaryExpr(smthng)).cwiseProduct(theta.unaryExpr(cosine)));
	res.col(1) = R.cwiseProduct((u.unaryExpr(smthng)).cwiseProduct(theta.unaryExpr(sine)));
	res.col(2) = R.cwiseProduct(u);

	X = X + res;
}