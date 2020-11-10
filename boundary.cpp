#include "boundary.hpp"
Boundary::Boundary(void){

    analytic = true;

}

void Boundary::_init_(pfbound fbound, pfstop fstop){

    Distance_Analytic = fbound;
    stop = fstop;
    analytic = true;
}

void Boundary::_init_(int dim, std::string lupbound, pfstop fstop){
    Distance_Numeric.Init(dim, lupbound);
    stop = fstop;
    analytic = false;
}

double Boundary::Dist(double* params, 
            Eigen::VectorXd & position, 
            Eigen::VectorXd & exitpoint,
            Eigen::VectorXd & normal){
    
    if(analytic){

        return Distance_Analytic(params,position,exitpoint,normal);
    }

    return Distance_Numeric.Eval(params, position, exitpoint, normal);
}
