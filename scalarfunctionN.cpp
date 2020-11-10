#include"scalarfunctionN.hpp"

ScalarFunctionN::ScalarFunctionN(void){
    function = Default_ScalarN;
    analytic = true;
}

void ScalarFunctionN::Init(pfscalarN input){
    function = input;
    analytic = true;
}

void ScalarFunctionN::Init(unsigned int dim, 
                          std::string input){
    lookuptable.Init(dim, input);
    analytic = false;
}

double ScalarFunctionN::Value(Eigen::VectorXd position, 
                            Eigen::VectorXd normal, 
                            double t){
    if(analytic){
        return function(position, normal, t);
    }

    return lookuptable.Eval(position);
}

double Default_ScalarN(Eigen::VectorXd X, Eigen::VectorXd N, double t){
    return 0.0f;
}
