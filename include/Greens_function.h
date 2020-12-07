#ifndef HARMONIC_GREENS_H
#define HARMONIC_GREENS_H

#define PI 3.1415926535897932384626433832795028884L

enum GREEN_FUNC_TYPE {
    HARMONIC,
    YUKAWA
};

#include <Eigen/Core>

void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G);

void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G,
        Eigen::VectorXd &G_int);

void Greens_function(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        const double &c,
        GREEN_FUNC_TYPE type,
        Eigen::VectorXd &G,
        Eigen::VectorXd &G_int,
        Eigen::MatrixXd &G_grad);

#endif //HARMONIC_GREENS_H
