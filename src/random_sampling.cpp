#include <random_sampling.h>

#include <igl/fast_winding_number.h>
#include <igl/slice_mask.h>

#include <random>
#include <chrono>

void random_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const int &n, Eigen::MatrixXd& P) {	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(-1., 1.);
	
	// Fill the bounding box with sampled points
	Eigen::MatrixXd Q = Eigen::MatrixXd::NullaryExpr(n, 3, [&]() { return dist(generator); });
	const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
	const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
	const Eigen::RowVector3d Vdiag = Vmax - Vmin;
	for (int q = 0; q < Q.rows(); q++) {
		Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
	}

	// using fast winding number to only get points in the mesh
	igl::FastWindingNumberBVH fwn_bvh;
	Eigen::VectorXf S;
	igl::fast_winding_number(V.cast<float>().eval(), F, 2, fwn_bvh);
	igl::fast_winding_number(fwn_bvh, 2, Q.cast<float>().eval(), S);
	igl::slice_mask(Q, S.array() > 0.5, 1, P);
}