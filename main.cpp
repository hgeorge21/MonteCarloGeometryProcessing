#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>

#include <timer.h>
#include <random_sampling.h>
#include <walk_on_spheres.h>

int main(int argc, char* argv[]) {
	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/bunny.off"), V, F);

	igl::opengl::glfw::Viewer viewer;
	viewer.data().point_size = 5;
	const int xid = viewer.selected_data_index;
	viewer.append_mesh();

	std::cout << R"(
    [space] Solves the PDE with current sampled points
    P       Show the sampled points
    p       Hide the sampled points
    +       Double the number of sample points
    -       Halve the number of sample points
    R,r     resets the sampled points
    G,g     changes the pre-listed boundary conditions
    )";
	std::cout << "\n";

    // sample the points
    int n_samples = 1048576;
    bool show_samples = true;
    Eigen::MatrixXd P;
    Eigen::VectorXd U;
    Eigen::VectorXd B;

	// TODO: add a few more boundary conditions
	int g_idx = 0;
	std::vector<const std::function<double(const Eigen::Vector3d&)>> boundary_funcs;
    boundary_funcs.emplace_back([](const Eigen::Vector3d &x) -> double { return x.maxCoeff(); });
    boundary_funcs.emplace_back([](const Eigen::Vector3d &x) -> double { return 2 * x.norm(); });

    std::vector<const std::string> boundary_strings;
    boundary_strings.emplace_back("||x||_1");
    boundary_strings.emplace_back("2 * ||x||_2");

    const auto &set_boundary = [&](const std::function<double(const Eigen::Vector3d&)> g) {
        B.resize(V.rows());
        for (int i = 0; i < V.rows(); i++) {
            B(i) = g(V.row(i));
        }
    };

	const auto update = [&]() {
	    Eigen::MatrixXd CM;
	    igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, V, -100, 100, CM);
	    viewer.data().set_colormap(CM);
	};

	const auto& set_points = [&]()
	{
		random_sampling(V, F, n_samples, P);
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
		if (show_samples)
			viewer.data_list[xid].set_points(P, (1. - (1. - orange.array()) * .8));
		else
			viewer.data_list[xid].clear_points();
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
		    std::cout << "Solving...\n";
			std::cout << "Solve time: " << timer(solve) << " sec\n";
			break;
		case '+':
			n_samples *= 2;
			set_points();
			std::cout << "# Sample points: " << P.rows() << "\n";
			break;
		case '-':
			n_samples /= 2;
			set_points();
			std::cout << "# Sample points: " << P.rows() << "\n";
			break;
		case 'P':
			show_samples = true;
			break;
		case 'p':
			show_samples = false;
			break;
        case 'G':
        case 'g':
            g_idx = (g_idx + 1) % boundary_funcs.size();
            set_boundary(boundary_funcs[g_idx]);
            std::cout << "Boundary condition: " << boundary_strings[g_idx] << "\n";
            break;
		case 'R':
		case 'r':
			set_points();
			break;
		default:
			return false;
		}
		update();
		return true;
	};

	viewer.data().set_mesh(V, F);
    set_boundary(boundary_funcs[g_idx]);
    std::cout << "Boundary condition: " << boundary_strings[g_idx] << "\n";
	std::cout << "Time used: " << timer(set_points) << " sec\n";

	update();
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}