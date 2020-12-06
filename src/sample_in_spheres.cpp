#include <sample_in_spheres.h>
#include <random>
#include <chrono>

#include <iostream>

void sample_in_spheres(const Eigen::MatrixXd& X, const Eigen::VectorXd& R, Eigen::MatrixXd& Y) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    std::uniform_real_distribution<double> dist_phi(0., 2 * PI);
    std::uniform_real_distribution<double> dist_costheta(-1., 1.);
    std::uniform_real_distribution<double> dist_u(0, 1);

    const auto& sine = [](double x) { return sin(x); };
    const auto& cosine = [](double x) { return cos(x); };

    Eigen::VectorXd phi   = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() {return dist_phi(generator); });
    Eigen::VectorXd theta = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return acos(dist_costheta(generator)); });
    Eigen::VectorXd r     = Eigen::VectorXd::NullaryExpr(X.rows(), [&]() { return std::cbrt(dist_u(generator)); });
    r = R.cwiseProduct(r);

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(X.rows(), 3);
    res.col(0) = r.cwiseProduct((theta.unaryExpr(sine)).cwiseProduct(phi.unaryExpr(cosine)));
    res.col(1) = r.cwiseProduct((theta.unaryExpr(sine)).cwiseProduct(phi.unaryExpr(sine)));
    res.col(2) = r.cwiseProduct(theta.unaryExpr(cosine));

    Y = X + res;
}