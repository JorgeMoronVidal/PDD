#include <eigen3/Eigen/Core>

double Equation_c(Eigen::VectorXd X, double t);
double Equation_f(Eigen::VectorXd X, double t);
double Equation_p(Eigen::VectorXd X, double t);
double Equation_g(Eigen::VectorXd X, double t);
double Equation_Psi(Eigen::VectorXd X, Eigen::VectorXd normal, double t);
double Equation_Varphi(Eigen::VectorXd X, Eigen::VectorXd normal, double t);
bool Stopping(Eigen::VectorXd X);
Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t);
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t);
Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t);
double Equation_u(Eigen::VectorXd X, double t);
