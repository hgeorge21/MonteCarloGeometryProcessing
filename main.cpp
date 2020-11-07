#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>

int main(int argc, char* argv[]) {
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/smthing.obj"), V, F);

	igl::opengl::glfw::Viewer viewer;
	std::cout << R"()";

	const auto update = [&]() {

	};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{

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