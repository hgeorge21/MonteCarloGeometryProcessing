#include <random_sampling.h>

#include <random>

#include <igl/AABB.h>
#include <igl/signed_distance.h>

void random_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& P) {
	Eigen::MatrixXd tmp;
	tmp.resize(10000, 3);
	
	Eigen::MatrixXd Q = Eigen::MatrixXd::Random(10000, 3);
	const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
	const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
	const Eigen::RowVector3d Vdiag = Vmax - Vmin;
	for (int q = 0; q < Q.rows(); q++)
	{
		Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
	}

	Eigen::VectorXd S;
	Eigen::VectorXi I;
	Eigen::MatrixXd C;
	Eigen::MatrixXd N;
	igl::signed_distance(Q, V, F, igl::SIGNED_DISTANCE_TYPE_DEFAULT, S, I, C, N);

	int count = 0;
	for (int i = 0; i < 10000; i++) {
		if (S(i) < 0) {
			tmp.row(count) = Q.row(i);
			count++;
		}
	}
	P = tmp.block(0, 0, 1000, 3);

}