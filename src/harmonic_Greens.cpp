#include <harmonic_Greens.h>

#include <iostream>

void harmonic_Greens(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& R,
                     Eigen::VectorXd &G, Eigen::VectorXd &G_int, Eigen::MatrixXd &G_grad) {
    Eigen::VectorXd r = (Y-X).rowwise().norm();
    G = 1./(4*PI) * (R - r).cwiseQuotient(R.cwiseProduct(r));
    G_int = R.cwiseProduct(R) / 6;

    const auto& inv_cube = [](double x) { return 1. / pow(x, 3); };
    Eigen::VectorXd tmp = (r.unaryExpr(inv_cube) - R.unaryExpr(inv_cube)) / (4 * PI);

    G_grad = Y - X;
    G_grad.col(0) = G_grad.col(0).cwiseProduct(tmp);
    G_grad.col(1) = G_grad.col(1).cwiseProduct(tmp);
    G_grad.col(2) = G_grad.col(2).cwiseProduct(tmp);
}