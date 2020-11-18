#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
#define kx 0.47
#define ky 0.89
#define C 2.0
double Equation_c(Eigen::VectorXd X, double t);
double Equation_f(Eigen::VectorXd X, double t);
double Equation_g(Eigen::VectorXd X, double t);
double Equation_p(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t);
Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_mu(Eigen::VectorXd X, double t);
double Equation_u(Eigen::VectorXd X, double t);
double Equation_Psi(Eigen::VectorXd X, Eigen::VectorXd normal, double t);
double Equation_Varphi(Eigen::VectorXd X, Eigen::VectorXd normal, double t);
double Equation_RBF(Eigen::VectorXd x , Eigen::VectorXd xj, double c2);
//bool Stopping(Eigen::VectorXd position);