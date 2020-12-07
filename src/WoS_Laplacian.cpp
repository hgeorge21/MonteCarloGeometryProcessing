#include <WoS_Laplacian.h>
#include <sample_on_spheres.h>
#include <interpolate.h>

#include <igl/AABB.h>

void WoS_Laplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P, Eigen::VectorXd& U) 
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const int max_itr = 5;
	const int n_walks = 64;

    Eigen::MatrixXd X;
    Eigen::MatrixXd C;
    Eigen::VectorXd R;
    Eigen::VectorXd D;
    Eigen::VectorXi I;
    U = Eigen::VectorXd::Zero(P.rows());

    const auto& handle_done_pts = [&](const Eigen::VectorXi &I_) {
        Eigen::Vector3d phi;
        for (int j = 0; j < I_.rows(); j++) {
            // interpolate between the three vertices and find value
            int f = I_(j);
            interpolate(C.row(j), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), phi);
            U(j) += phi(0) * B(F(f,0)) + phi(1) * B(F(f,1)) + phi(2) * B(F(f,2));
        }
    };

	// start the random walks
	for (int i = 0; i < n_walks; i++) {
		X = P;
        int itr = 0;
		while (itr < max_itr) {
			aabb.squared_distance(V, F, X, D, I, C);
            R = D.cwiseSqrt();
			sample_on_spheres(R, X);
			itr++;
		}
		handle_done_pts(I);
	}
	U = U / n_walks;
}