#include <WoS_Poisson.h>
#include <sample_on_spheres.h>
#include <sample_in_spheres.h>
#include <interpolate.h>
#include <Greens_function.h>

#include <igl/AABB.h>

void WoS_Poisson(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
                 std::function<double(const Eigen::Vector3d &)> f,
                 const bool &use_pt_src, const Eigen::RowVector3d& point_source,
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
    // ================

    // start the random walk
    for (int i = 0; i < n_walks; i++) {
        X = P;
        int itr = 0;
        while (itr < max_itr) {
            // find closest point on the boundary and change radius
            aabb.squared_distance(V, F, X, D, I, C);
            R = D.cwiseSqrt();

            if(use_pt_src)
                sample_in_spheres(X, R, point_source, Y);
            else
                sample_in_spheres(X, R, Y);
            vols = R.unaryExpr([](double r)->double{ return 4.*PI/3*pow(r, 3); });
            Greens_function(X, Y, R, 0, HARMONIC, G);

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