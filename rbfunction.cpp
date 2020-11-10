#include "rbfunction.hpp"

RBFunction::RBFunction(void){
    function = Default_RBF;
}

void RBFunction::Init(pRBF input){
    function = input;
}


double RBFunction::Value (Eigen::VectorXd x, 
                          Eigen::VectorXd x_i,
                          double c2){
    return function(x, x_i, c2);
}

double Default_RBF(Eigen::VectorXd x, 
                  Eigen::VectorXd x_i,
                  double c2){
    return 0.0f;
}