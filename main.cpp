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
	const int xid = viewer.selected_data_index;
    viewer.data().point_size = 5;
	viewer.append_mesh();
	const int yid = viewer.selected_data_index;
    viewer.data().point_size = 8;

	std::cout << R"(
    [space] Solves the PDE with current sampled points
    H,h     Hide boundary points
    S,s     Show boundary points
    +       Double the number of sample points
    -       Halve the number of sample points
    R,r     resets the sampled points
    G,g     changes the pre-listed boundary conditions
    )";
	std::cout << "\n";

    // sample the points
    int n_samples = 1048576;
    bool show_boundary = true;
    double solve_time;
    Eigen::MatrixXd P;
    Eigen::VectorXd U;
    Eigen::VectorXd B;
    Eigen::MatrixXd BCM;

	// TODO: add a few more boundary conditions
	int g_idx = 0;
	std::vector<std::function<double(const Eigen::Vector3d&)>> boundary_funcs;
    boundary_funcs.emplace_back([](const Eigen::Vector3d &x) -> double { return x(0); });
    boundary_funcs.emplace_back([](const Eigen::Vector3d &x) -> double { return 2 * x.norm(); });
    boundary_funcs.emplace_back([](const Eigen::Vector3d &x) -> double { return 5.; });

    std::vector<std::string> boundary_strings;
    boundary_strings.emplace_back("g(x) = x_0");
    boundary_strings.emplace_back("g(x) = 2 * ||x||");
    boundary_strings.emplace_back("g(x) = 5");

    const auto &set_boundary = [&](std::function<double(const Eigen::Vector3d&)> g) {
        B.resize(V.rows());
        for (int i = 0; i < V.rows(); i++) {
            B(i) = g(V.row(i));
        }
        igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, V, B.minCoeff(), B.maxCoeff(), BCM);
        viewer.data_list[yid].set_points(V, BCM);
    };

	const auto& set_points = [&]()
	{
		random_sampling(V, F, n_samples, P);
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
        viewer.data_list[xid].set_points(P, (1. - (1. - orange.array()) * .8));
	};
	const auto& show_bndry_pts = [&]()
    {
        if (show_boundary)
            viewer.data_list[yid].set_points(V, BCM);
        else
            viewer.data_list[yid].clear_points();
    };

	const auto& solve = [&]()
	{
		walk_on_spheres(V, F, B, P, U);
		Eigen::MatrixXd CM;
		igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, U, U.minCoeff(), U.maxCoeff(), CM);
		viewer.data_list[xid].set_points(P, CM);
	};

    const auto update = [&]() {
        show_bndry_pts();
    };

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
		    std::cout << "Solving...\n";
		    solve_time = timer(solve);
			std::cout << "Solve time: " << solve_time << " sec\n";
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
        case 'S':
		case 's':
			show_boundary = true;
			break;
		case 'H':
		case 'h':
			show_boundary = false;
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
	set_points();
    std::cout << "# Sample points: " << P.rows() << "\n";

	update();
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}