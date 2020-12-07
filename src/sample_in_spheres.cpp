#include <sample_in_spheres.h>
#include <random>
#include <chrono>

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

    // This sampling ensures uniform distribution within sphere
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(X.rows(), 3);
    res.col(0) = r.cwiseProduct((theta.unaryExpr(sine)).cwiseProduct(phi.unaryExpr(cosine)));
    res.col(1) = r.cwiseProduct((theta.unaryExpr(sine)).cwiseProduct(phi.unaryExpr(sine)));
    res.col(2) = r.cwiseProduct(theta.unaryExpr(cosine));

    Y = X + res;
}


void sample_in_spheres(const Eigen::MatrixXd& X, const Eigen::VectorXd& R, const Eigen::RowVector3d &point_source,
                       Eigen::MatrixXd& Y) {
    sample_in_spheres(X, R, Y);
    for(int i = 0; i < X.rows(); i++) {
        if((point_source - X.row(i)).norm() < R(i))
            Y.row(i) = point_source;
    }
}

