#include <walk_on_spheres.h>
#include <pts_on_spheres.h>
#include <interpolate.h>

#include <igl/AABB.h>

void walk_on_spheres(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P, Eigen::VectorXd& U) 
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const float eps = 1e-6;
	const int max_itr = 20;
	const int n_walks = 120;

	U = Eigen::VectorXd::Zero(P.rows());
	
	// start the random walk
	for (int i = 0; i < n_walks; i++) {
		Eigen::MatrixXd X = P;
		Eigen::VectorXd R = std::numeric_limits<double>::max() * Eigen::VectorXd::Ones(P.rows());

		int itr = 0;
		Eigen::VectorXd sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		while (itr < max_itr) {
			// find closest point on the boundary and change radius
			aabb.squared_distance(V, F, P, sqrD, I, C);

			R = R.cwiseMin(sqrD);
			pts_on_spheres(R, X);
			itr++;
		}
		for (int j = 0; j < I.rows(); j++) {
			// interpolate between the three vertices and find value
			Eigen::Vector3d phi = Eigen::Vector3d::Ones() / 3;
			//interpolate(C.row(j), V.row(F(I(j), 0)), V.row(F(I(j), 1)), V.row(F(I(j), 2)), phi);

			U(j) += phi(0) * B(I(j),0) + phi(1) * B(I(j),1) + phi(2) * B(I(j),2);
		}
	}
	U = U / n_walks;
}