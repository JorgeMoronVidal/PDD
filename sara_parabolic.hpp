#include <eigen3/Eigen/Core>

float Equation_c(Eigen::VectorXf X, float t);
float Equation_f(Eigen::VectorXf X, float t);
float Equation_g(Eigen::VectorXf X, float t);
float Equation_p(Eigen::VectorXf X, float t);
Eigen::VectorXf Equation_b(Eigen::VectorXf X, float t);
Eigen::VectorXf Equation_F(Eigen::VectorXf X, float t);
Eigen::VectorXf Equation_mu(Eigen::VectorXf X, float t);
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, float t);
float Equation_u(Eigen::VectorXf X, float t);
float Equation_RBF(Eigen::VectorXf X , Eigen::VectorXf Xj, float c2);
bool Stopping(Eigen::VectorXf position);