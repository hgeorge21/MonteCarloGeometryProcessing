#ifndef HARMONIC_GREENS_H
#define HARMONIC_GREENS_H

#define PI 3.1415926535897932384626433832795028884L

#include <Eigen/Core>

void harmonic_Greens(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Y,
        const Eigen::VectorXd& R,
        Eigen::VectorXd &G,
        Eigen::VectorXd &G_int,
        Eigen::MatrixXd &G_grad);

#endif //HARMONIC_GREENS_H
