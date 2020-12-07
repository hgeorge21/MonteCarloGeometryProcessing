#ifndef BOUNDARY_SETUP_H
#define BOUNDARY_SETUP_H

#include <Eigen/Core>
#include <WoS_biharmonic.h>
#include <WoS_Laplacian.h>
#include <WoS_Poisson.h>

typedef std::pair<std::string, std::function<double(const Eigen::Vector3d&)>> f_pairs;

void boundary_setup(
        const Eigen::Vector3d &source_point,
        std::vector<std::pair<std::string, std::vector<f_pairs>>> &solver_funcs) {

    // Laplacian PDE boundary functions
    std::vector<f_pairs> Laplacian_boundary_funcs;
    Laplacian_boundary_funcs.emplace_back(
            std::make_pair("g(x) = 1 / ||x||",
                           [](const Eigen::Vector3d &x) -> double { return 1 / x.norm(); }));
    Laplacian_boundary_funcs.emplace_back(
            std::make_pair("g(x) = x_0",
                           [](const Eigen::Vector3d &x) -> double { return x(0); }));
    Laplacian_boundary_funcs.emplace_back(
            std::make_pair("g(x) = 2 * ||x||",
                           [](const Eigen::Vector3d &x) -> double { return 2 * x.norm(); }));

    // Poisson PDE boundary functions
    std::vector<f_pairs> Poisson_boundary_funcs;
    Poisson_boundary_funcs.emplace_back(
            std::make_pair("g(x) = 1 / ||x - (source_point)||",
                           [&](const Eigen::Vector3d &x) -> double { return 1 / (x - source_point).norm(); }));
    
    // two boundary functions for biharmonic equation
    std::vector<f_pairs> biharmonic_boundary_funcs;
    biharmonic_boundary_funcs.emplace_back(
            std::make_pair("g(x) = 1 / ||x - (source_point)||",
                           [&](const Eigen::Vector3d &x) -> double { return 1 / (x - source_point).norm(); }));
    biharmonic_boundary_funcs.emplace_back(
            std::make_pair("h(x) = x_0",
                           [&](const Eigen::Vector3d &x) -> double { return x(0); }));

    // add the function handlers
    solver_funcs.emplace_back(
            std::make_pair("Laplacian", Laplacian_boundary_funcs));
    solver_funcs.emplace_back(
            std::make_pair("Poisson", Poisson_boundary_funcs));
    solver_funcs.emplace_back(
            std::make_pair("biharmonic", biharmonic_boundary_funcs));
}

#endif //BOUNDARY_SETUP_H
