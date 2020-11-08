#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/parula.h>
#include <igl/colormap.h>

#include <iostream>

int main(int argc, char* argv[]) {
	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/smthing.obj"), V, F);

	// set up AABB tree
	igl::AABB<Eigen::MatrixXd, 3> AABB_tree;
	AABB_tree.init(V, F);

	igl::opengl::glfw::Viewer viewer;
	std::cout << R"()";

	double scale = 100.;
	bool show_org_mesh = true;

	const auto update = [&]() {
		if (show_org_mesh) {
			Eigen::MatrixXd C;
			igl::parula(V, -scale, scale, C);
			viewer.data().set_colors(C);
		}
		else {
			Eigen::MatrixXd CM;
			igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, V, -scale, scale, CM);
			viewer.data().set_colormap(CM);
		}
	};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
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

	update();
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}