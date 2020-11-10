#include "scalarfunction.hpp"

ScalarFunction::ScalarFunction(void){
    function = Default_Scalar;
    analytic = true;
}

void ScalarFunction::Init(pfscalar input){
    function = input;
    analytic = true;
}

void ScalarFunction::Init(int dim, 
                          std::string input){
    lookuptable.Init(dim, input);
    analytic = false;
}

double ScalarFunction::Value(Eigen::VectorXd position, 
                            double t){
    if(analytic){
        return function(position, t);
    }

    return lookuptable.Eval(position);
}

double Default_Scalar(Eigen::VectorXd X, double t){
    return 0.0;
}
