#include <walk_on_sphere.h>
#include <random>
#include <functional>

#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>


void walk_on_sphere(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const igl::AABB<Eigen::MatrixXd, 3> aabb,
	const std::function<float(Eigen::Vector3d, double, bool)> u_hat, const Eigen::Vector3d &x0, double &u) 
{
	const float eps = 1e-6;
	const int max_itr = 40;
	const int n_walks = 120;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> dist_theta(0., 2 * std::_Pi);
	std::uniform_real_distribution<double> dist_phi(0., 2 * std::_Pi);
	
	auto& get_sphere_vec = [&](const double& theta, const double& phi) {
		Eigen::Vector3d vec;
		vec << sin(theta) * cos(phi),
			sin(theta)* sin(phi),
			cos(theta);
		return vec;
	};

	u = 0;
	for (int i = 0; i < n_walks; i++) {
		double r = std::numeric_limits<double>::max();
		Eigen::RowVector3d x = x0;

		int itr = 0;
		while (itr < max_itr && r > eps) {
			// find closest point on the boundary and change radius
			int idx;
			Eigen::RowVector3d c;
			double d = aabb.squared_distance(V, F, x, d, idx, c);
			r = std::min(r, std::sqrt(d));
			
			// add to u if needed
			u += u_hat(x, r, false);
			
			x = x + r * get_sphere_vec(dist_theta(generator), dist_phi(generator));
			itr++;
		}
		u += u_hat(x, r, true);
	}
	u = u / n_walks;
}