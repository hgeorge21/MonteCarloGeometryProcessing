#include <random_sampling.h>

#include <igl/fast_winding_number.h>
#include <igl/slice_mask.h>

void random_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const int &n, Eigen::MatrixXd& P) {
	Eigen::MatrixXd tmp;
	tmp.resize(n, 3);
	
	Eigen::MatrixXd Q = Eigen::MatrixXd::Random(n, 3);
	const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
	const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
	const Eigen::RowVector3d Vdiag = Vmax - Vmin;
	for (int q = 0; q < Q.rows(); q++) {
		Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
	}

	igl::FastWindingNumberBVH fwn_bvh;
	Eigen::VectorXf S;
	igl::fast_winding_number(V.cast<float>().eval(), F, 2, fwn_bvh);
	igl::fast_winding_number(fwn_bvh, 2, Q.cast<float>().eval(), S);
	igl::slice_mask(Q, S.array() > 0.5, 1, P);
}