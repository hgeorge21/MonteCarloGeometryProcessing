#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/colormap.h>

#include <iostream>

#include <random_sampling.h>
#include <walk_on_spheres.h>

int main(int argc, char* argv[]) {
	// timer function 
	// (copied from https://github.com/libigl/libigl/blob/master/tutorial/717_FastWindingNumber/main.cpp?fbclid=IwAR2sUXck8cvuTrhakA2NGNCvDSXxk04sTMJtUvLtGrZc3vFl_dRJ9TGkn3k)
	const auto time = [](std::function<void(void)> func)->double
	{
		const double t_before = igl::get_seconds();
		func();
		const double t_after = igl::get_seconds();
		return t_after - t_before;
	};


	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::read_triangle_mesh((argc > 1 ? argv[1] : "../../../data/bunny.off"), V, F);

	// compute the boundary conditions
	// for simplicity, use 2||x|| for now
	Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		B(i) = 2 * V.row(i).norm();
	}

	igl::opengl::glfw::Viewer viewer;
	viewer.data().point_size = 5;
	const int xid = viewer.selected_data_index;
	viewer.append_mesh();
	const int yid = viewer.selected_data_index;

	std::cout << R"()";

	double scale = 100.;
	bool show_org_mesh = true;

	const auto update = [&]() {};

	// sample the points
	int n_samples = 1048576;
	bool show_samples = true;
	Eigen::MatrixXd P;
	Eigen::VectorXd U;
	
	const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
	const auto& set_points = [&]()
	{
		random_sampling(V, F, n_samples, P);
		if (show_samples)
		{
			viewer.data_list[xid].set_points(P, (1. - (1. - orange.array()) * .8));
		}
		else
		{
			viewer.data_list[xid].clear_points();
		}
	};

	const auto& solve = [&]()
	{
		walk_on_spheres(V, F, B, P, U);
		Eigen::MatrixXd CM;
		igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, U, U.minCoeff(), U.maxCoeff(), CM);
		viewer.data_list[xid].set_points(P, CM);
	};


	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
			std::cout << time(solve) << "\n";
			break;
		case '+':
			n_samples *= 2;
			set_points();
			std::cout << "Samples used: " << P.rows() << "\n";
			break;
		case '-':
			n_samples /= 2;
			set_points();
			std::cout << "Samples used: " << P.rows() << "\n";
		case 'P':
			show_samples = true;
			break;
		case 'p':
			show_samples = false;
			break;
		case 'R':
		case 'r':
			set_points();
			break;
		case 'M':
		case 'm':
			show_org_mesh ^= 1;
			break;
		default:
			return false;
		}
		update();
		return true;
	};

	viewer.data().set_mesh(V, F);
	std::cout << "Time used: " << time(set_points) << "s\n";

	update();
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}