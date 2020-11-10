#include "sara.hpp"
#include <math.h>
#include <iostream>

double   Equation_c(Eigen::VectorXd X, double t){
    return -(pow(X(0),2.0)+pow(X(1),2.0))*exp(-X(0));
}
double   Equation_f(Eigen::VectorXd X, double t){
    double aux = (pow(X(0),2.0)+pow(X(1),2.0))*exp(-X(0));
    //In order to enforce the solution to be positive
    //aa = 1.5f, bb =2.0f, C = 2.1f
    //res  = res= (aa^2+aux).*cos(aa*x) + (bb^2+aux).*sin(bb*y) + aa*sin(aa*x) + C*aux;
    return (2.25f + aux)*cos(1.5f*X(0)) + (4.0f + aux)*
    sin(2.0f*X(1)) + 1.5f*sin(1.5f*X(0)) + 2.1f*aux;
}
  double Equation_g(Eigen::VectorXd X, double t){
    return Equation_u(X,t);
}
  Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t){
    Eigen::VectorXd b;
    b = Eigen::VectorXd::Constant(2,1.0f);
    b(1) = 0.0f;
    return b;
}
  Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t){
    Eigen::VectorXd F(2);
    F = Eigen::VectorXd::Zero(2);
    return F;
}
  Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t){
    Eigen::MatrixXd sigma(2,2);
    sigma << 1.41421356237f ,0.0f,
            0.0f, 1.41421356237f;
    return sigma;
}
  double Equation_u(Eigen::VectorXd X, double t){
    // res= cos(aa*x) + sin(bb*y) + C;
    return cos(1.5f * X(0)) + sin(2.0f*X(1)) + 2.1f;
}

double Equation_RBF(Eigen::VectorXd X , Eigen::VectorXd Xj, double c2){
    double r2 = pow(X(0)-Xj(0),2) + pow(X(1)-Xj(1),2);
    return sqrt(r2 + c2);
}
bool Stopping(Eigen::VectorXd position){
  return true;
}