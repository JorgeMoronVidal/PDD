#include "edepaco.hpp"

double Equation_c(Eigen::VectorXd X, double t){
    return -pow(X.norm(),2)* exp(-X(0));
}
double Equation_f(Eigen::VectorXd X, double t){
    return 2.0 * sin(5.0*t + 7.0*X(0) + 9.0*X(1)) + 130.0*
    cos(5.0*t + 7.0*X(0) + 9.0*X(1)) -Equation_c(X,t) * Equation_u(X,t);
}
double Equation_p(Eigen::VectorXd X, double t){
    return 2.1f + cos(7.0 * X(0) + 9.0 * X(1));
}
double Equation_g(Eigen::VectorXd X, double t){
    return 2.1f + cos(5.0 * t + 7.0 * X(0) + 9.0 * X(1));
}
double Equation_Psi(Eigen::VectorXd X, Eigen::VectorXd normal, double t){
    return -(normal(0) * 7.0 + normal(1) * 9.0) * sin(5.0 * t + 7.0 * X(0) + 9.0 * X(1));
            //+ Equation_u(X,t);
    //return 0.0;
}
double Equation_Varphi(Eigen::VectorXd X, Eigen::VectorXd normal, double t){
    return -0.0;
    //return (-(normal(0) * 7.0 + normal(1) * 9.0) * sin(5.0 * t + 7.0 * X(0) + 9.0 * X(1)))/Equation_u(X,t);
}
bool Stopping(Eigen::VectorXd X){
    return false;
}
Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t){
    Eigen::VectorXd b(X.size());
    b << 1.0, 0.0;
    return b;
}
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t){
    Eigen::VectorXd F(X.size());
    F << 7.0, 9.0;
    F = -F*sin(5.0 *t + 7.0 * X(0) + 9.0 * X(1) );
    F = -Equation_sigma(X,t).transpose()*F;
    return F;
}
Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t){
    Eigen::MatrixXd sigma(X.size(), X.size());
    sigma << 1.41421356237, 0.0,
             0.0, 1.41421356237;
    return sigma;
}
double Equation_u(Eigen::VectorXd X, double t){
    return 2.1 + cos(5.0 * t + 7.0 * X(0) +
    9.0 * X(1));
}
