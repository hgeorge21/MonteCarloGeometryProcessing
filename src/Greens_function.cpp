#include <Greens_function.h>


void Greens_function(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& R,
                     const double &c, GREEN_FUNC_TYPE type, Eigen::VectorXd &G) {
    Eigen::VectorXd r = (Y - X).rowwise().norm();

    if(type == HARMONIC) {
        G = 1. / (4 * PI) * (R - r).cwiseQuotient(R.cwiseProduct(r));
    } else if(type == YUKAWA) {
        Eigen::VectorXd tmp = sqrt(c) * (R-r);
        G = 1. / (4 * PI) * tmp.unaryExpr([](double x){ return sinh(x); });
        G = G.cwiseQuotient(r.cwiseProduct((sqrt(c)*R).unaryExpr([](double x){ return sinh(x); })));
    }
}


void Greens_function(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& R,
                     const double &c, GREEN_FUNC_TYPE type, Eigen::VectorXd &G, Eigen::VectorXd &G_int) {
    Greens_function(X, Y, R, c, type, G);
    if(type == HARMONIC) {
        G_int = R.cwiseProduct(R) / 6;
    } else if(type == YUKAWA){
        G_int = R.unaryExpr([&](double r) {
            return (1-sqrt(c)*r/sinh(sqrt(c)*r)) / c;
        });
    }
}


void Greens_function(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& R,
                     const double &c, GREEN_FUNC_TYPE type,
                     Eigen::VectorXd &G, Eigen::VectorXd &G_int, Eigen::MatrixXd &G_grad) {
    Eigen::VectorXd r = (Y - X).rowwise().norm();
    if(type == HARMONIC) {
        Greens_function(X, Y, R, c, type, G, G_int);

        const auto &inv_cube = [](double x) { return 1. / pow(x, 3); };
        Eigen::VectorXd tmp = (r.unaryExpr(inv_cube) - R.unaryExpr(inv_cube)) / (4 * PI);
        G_grad = Y - X;
        G_grad.col(0) = G_grad.col(0).cwiseProduct(tmp);
        G_grad.col(1) = G_grad.col(1).cwiseProduct(tmp);
        G_grad.col(2) = G_grad.col(2).cwiseProduct(tmp);
    } else if(type == YUKAWA) {
        const auto& sinh_f = [](double x){ return sinh(x); };
        const auto& cosh_f = [](double x){ return cosh(x); };

        G_grad = (Y-X) / (4*PI);
        Eigen::VectorXd tmp1, tmp2;
        tmp1 = (r.cwiseInverse() - R.cwiseInverse()) * sqrt(c) * (sqrt(c)*(R-r)).unaryExpr(cosh_f).cwiseQuotient(
                r.cwiseProduct((sqrt(c)*R).unaryExpr(sinh_f)));
        tmp2 = ((sqrt(c)*(R-r)).unaryExpr(sinh_f)).cwiseQuotient(r.cwiseProduct((sqrt(c)*R).unaryExpr(sinh_f))).cwiseProduct(
                r.unaryExpr([](double x){return 1./(x*x); }) + sqrt(c)*((sqrt(c)*R.unaryExpr(cosh_f))).cwiseQuotient(
                        R.unaryExpr([&](double x) { return x*sinh(x*sqrt(c)); })));
        tmp1 = tmp1 + tmp2;
        G_grad.col(0) = G_grad.col(0).cwiseProduct(tmp1);
        G_grad.col(1) = G_grad.col(1).cwiseProduct(tmp1);
        G_grad.col(2) = G_grad.col(2).cwiseProduct(tmp1);
    }
}