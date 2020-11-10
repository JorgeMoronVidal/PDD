#include <eigen3/Eigen/Core>

double Equation_c(Eigen::VectorXd X, double t);
double Equation_f(Eigen::VectorXd X, double t);
double Equation_g(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t);
Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t);
double Equation_u(Eigen::VectorXd X, double t);
double Equation_RBF(Eigen::VectorXd X , Eigen::VectorXd Xj, double c2);
bool Stopping(Eigen::VectorXd position);