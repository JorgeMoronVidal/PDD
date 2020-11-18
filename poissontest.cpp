#include "poissontest.hpp"
double Equation_c(Eigen::VectorXd X, double t){
    return 0.0;
}
double Equation_f(Eigen::VectorXd X, double t){
    return M_PI*M_PI*(kx*kx + ky*ky)*sin(kx*M_PI*X(0))*sin(ky*M_PI*X(1));
}
double Equation_g(Eigen::VectorXd X, double t){
    return Equation_u(X,t);
}
double Equation_p(Eigen::VectorXd X, double t){
    return Equation_u(X,t);
}
Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t){
    return 0.0 * X;
}
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t){
    Eigen::VectorXd F;
    F.resize(X.size());
    F(0) =  kx*M_PI*cos(kx*M_PI*X(0))*sin(ky*M_PI*X(1));
    F(1) =  ky*M_PI*sin(kx*M_PI*X(0))*cos(ky*M_PI*X(1));
    F = -Equation_sigma(X,t).transpose()*F;
    return F;
}
Eigen::VectorXd Equation_mu(Eigen::VectorXd X, double t){
    return Equation_F(X,t)/Equation_u(X,t);
}
Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t){
    /*Eigen::MatrixXd sigma;
    sigma.resize(2,2);
    sigma(0,0) = 1.4142136;
    sigma(1,1) = 1.4142136;
    return sigma;
    */
    return Eigen::MatrixXd::Identity(2,2)*1.4142136;
}
double Equation_Psi(Eigen::VectorXd X, Eigen::VectorXd normal, double t){
    return normal(0)*kx*M_PI*cos(kx*M_PI*X(0))*sin(ky*M_PI*X(1))+
           normal(1)*ky*M_PI*sin(kx*M_PI*X(0))*cos(ky*M_PI*X(1))+
           Equation_u(X,t);
}
double Equation_Varphi(Eigen::VectorXd X, Eigen::VectorXd normal, double t){
    return -1;
}
double Equation_u(Eigen::VectorXd X, double t){
    return sin(kx*M_PI*X(0))*sin(ky*M_PI*X(1)) + C;
}
double Equation_RBF(Eigen::VectorXd x , Eigen::VectorXd xj, double c2){
    double r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return sqrt(r2 + c2);
}
/*bool Stopping(Eigen::VectorXd position){
    bool stop = true;
    return stop;
}*/