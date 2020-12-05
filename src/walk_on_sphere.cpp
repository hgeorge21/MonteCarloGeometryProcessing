#include <walk_on_sphere.h>
#include <random>
#include <functional>

#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>


void walk_on_spheres(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P, Eigen::VectorXd& U) 
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const float eps = 1e-6;
	const int max_itr = 40;
	const int n_walks = 120;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> dist_theta(0., 2 * std::_Pi);
	std::uniform_real_distribution<double> dist_phi(0., 2 * std::_Pi);
	
	// get point on a sphere
	auto& get_pts_on_sphere = [&](const Eigen::VectorXd &R) {
		Eigen::MatrixXd Vec = Eigen::MatrixXd::Zero(P.rows(), 3);

		for (int k = 0; k < P.rows(); k++) {
			double theta = dist_theta(generator);
			double phi = dist_phi(generator);

			Vec.row(k) << sin(theta) * cos(phi), sin(theta)* sin(phi), cos(theta);
			Vec.row(k) *= R(k);
		}
		return Vec;
	};

	U = Eigen::VectorXd::Zero(P.rows());
	Eigen::MatrixXd X = P;

	// start the random walk
	for (int i = 0; i < n_walks; i++) {
		Eigen::VectorXd R = std::numeric_limits<double>::max() * Eigen::VectorXd::Ones(P.rows());

		int itr = 0;
		Eigen::VectorXd sqrD;
		Eigen::VectorXd I;
		Eigen::MatrixXd C;
		while (itr < max_itr) {
			// find closest point on the boundary and change radius
			aabb.squared_distance(V, F, P, sqrD, I, C);
			
			R = R.cwiseMin(sqrD);
			X = X + get_pts_on_sphere(R);
			itr++;
		}
		for (int j = 0; j < I.rows(); j++) {
			// interpolate between the three vertices and find value
			// TODO: use barycentric for now...
			U(j) += (B(I(j),0)+B(I(j),1)+B(I(j),2)) / 3;
		}
	}
	U = U / n_walks;
}