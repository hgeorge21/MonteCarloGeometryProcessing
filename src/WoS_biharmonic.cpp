#include <WoS_biharmonic.h>
#include <sample_on_spheres.h>
#include <sample_in_spheres.h>
#include <Greens_function.h>
#include <interpolate.h>

#include <igl/AABB.h>

void WoS_biharmonic(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                    const Eigen::VectorXd& Bg, const Eigen::VectorXd& Bh,
                    const bool &use_pt_src, const Eigen::RowVector3d& point_source,
                    const Eigen::MatrixXd& P, Eigen::VectorXd& U) {
    igl::AABB<Eigen::MatrixXd, 3> aabb;
    aabb.init(V, F);

    const int max_itr = 5;
    const int n_walks = 64;

    Eigen::VectorXd W;
    Eigen::MatrixXd X;
    Eigen::MatrixXd Y;
    Eigen::VectorXd R;

    Eigen::MatrixXd CX;
    Eigen::VectorXd DX;
    Eigen::VectorXi IX;
    Eigen::MatrixXd CY;
    Eigen::VectorXd DY;
    Eigen::VectorXi IY;

    U = Eigen::VectorXd::Zero(P.rows());
    W = Eigen::VectorXd::Zero(P.rows());

    const auto &handle_boundary_pts = [&](const Eigen::MatrixXd &B_, const Eigen::MatrixXd &C_,
            const Eigen::VectorXi &I_, Eigen::VectorXd &U_) {
        Eigen::Vector3d phi;
        for (int j = 0; j < I_.rows(); j++) {
            // interpolate between the three vertices and find value
            int f = I_(j);
            interpolate(C_.row(j), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), phi);
            U_(j) += phi(0) * B_(F(f, 0)) + phi(1) * B_(F(f, 1)) + phi(2) * B_(F(f, 2));
        }
    };

    // additional parts
    Eigen::VectorXd vols = Eigen::VectorXd::Zero(P.rows());
    Eigen::VectorXd res = Eigen::VectorXd::Zero(P.rows());
    Eigen::VectorXd G;
    // ================

    // start the random walk
    for (int i = 0; i < n_walks; i++) {
        X = P;
        int itr = 0;
        while (itr < max_itr) {
            // find closest point on the boundary and change radius
            aabb.squared_distance(V, F, X, DX, IX, CX);
            R = DX.cwiseSqrt();

            // walk y once
            if(use_pt_src)
                sample_in_spheres(X, R, point_source, Y);
            else
                sample_in_spheres(X, R, Y);
            aabb.squared_distance(V, F, Y, DY, IY, CY);
            handle_boundary_pts(Bh, CY, IY, W);

            vols = R.unaryExpr([](double r)->double{ return 4.*PI/3*pow(r, 3); });
            Greens_function(X, Y, R, 0, HARMONIC, G);

            U = U + vols.cwiseProduct(W).cwiseProduct(G);
            sample_on_spheres(R, X);
            itr++;
        }
        // handle the rest of the points
        handle_boundary_pts(Bg, CX, IX, U);
    }
    U = U / n_walks;
}