#include <WoS_Poisson.h>
#include <sample_on_spheres.h>
#include <sample_in_spheres.h>
#include <interpolate.h>
#include <harmonic_Greens.h>

#include <igl/AABB.h>
#include <numeric>

void WoS_Poisson(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
                 std::function<double(const Eigen::Vector3d &)> f,
                 const Eigen::MatrixXd& P, Eigen::VectorXd& U) {
    igl::AABB<Eigen::MatrixXd, 3> aabb;
    aabb.init(V, F);

    const int max_itr = 5;
    const int n_walks = 64;

    Eigen::MatrixXd X;
    Eigen::MatrixXd C;
    Eigen::VectorXd D;
    Eigen::VectorXd R;
    Eigen::VectorXi I;
    U = Eigen::VectorXd::Zero(P.rows());

    const auto &handle_pts = [&](const Eigen::VectorXi &I_) {
        Eigen::Vector3d phi;
        for (int j = 0; j < I_.rows(); j++) {
            // interpolate between the three vertices and find value
            int f = I_(j);
            interpolate(C.row(j), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), phi);
            U(j) += phi(0) * B(F(f, 0)) + phi(1) * B(F(f, 1)) + phi(2) * B(F(f, 2));
        }
    };

    // additional parts
    Eigen::VectorXd vols = Eigen::VectorXd::Zero(P.rows());
    Eigen::VectorXd Y_res = Eigen::VectorXd::Zero(P.rows());
    Eigen::MatrixXd Y;
    Eigen::VectorXd G;
    Eigen::VectorXd G_int;
    Eigen::MatrixXd G_grad;
    // ================

    // start the random walk
    for (int i = 0; i < n_walks; i++) {
        X = P;
        R = std::numeric_limits<double>::max() * Eigen::VectorXd::Ones(P.rows());

        int itr = 0;
        while (itr < max_itr) {
            // find closest point on the boundary and change radius
            aabb.squared_distance(V, F, X, D, I, C);
            R = D.cwiseSqrt();

            sample_in_spheres(X, R, Y);
            vols = R.unaryExpr([](double x)->double{ return 4.*PI/3*pow(x, 3); });
            harmonic_Greens(X, Y, R, G, G_int, G_grad);

            for(int k = 0; k < Y.rows(); k++)
                Y_res(k) = f(Y.row(k));
            U = U + vols.cwiseProduct(Y_res).cwiseProduct(G);
            sample_on_spheres(R, X);
            itr++;
        }
        // handle the rest of the points
        handle_pts(I);
    }
    U = U / n_walks;
}
//
//double k = std::cbrt(uniform01(generator));
//double theta_y = 2 * 3.14 * uniform01(generator);
//double phi_y = 3.14 * uniform01(generator);
//double xx = sin(phi_y) * cos(theta_y);
//double yy = sin(phi_y) * cos(theta_y);
//double zz = cos(phi_y);
//
//Eigen::RowVector3d sample_y(xx, yy, zz);
//sample_y = sample * k * radius + center;
//
//// double fy = 1;
//double r2 = (sample_y - sourcePoint).squaredNorm();
//double fy = c * std::pow(exp(1.0), -r2);
//
//double volume = 1.0/3.0 * radius * radius; //  4 * pi *radius canceled by G(x, y)
//// (R-r)/rR = (1-k)R/kRR = (1-k)/(k * radius)
//U(i) += volume * fy * (1-k) / k;