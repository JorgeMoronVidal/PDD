#include"vectorfunction.hpp"

VectorFunction::VectorFunction(void){
    function = Default_Vector;
    analytic = true;
}

void VectorFunction::Init(pfvector input){
    function = input;
    analytic = true;
}

void VectorFunction::Init(int dim, 
                          std::string input){
    /*std::string aux;
    lookuptable.resize(dim);
    for(int i = 0; i < dim; i++){
        aux = input + "/" + std::to_string(i) +"_component";
        lookuptable[i].Init(dim, aux);
        aux.clear();
    }
    analytic = false;
    */
   lookuptable.Init(dim, input);
   analytic = false;
}

Eigen::VectorXd VectorFunction::Value(Eigen::VectorXd position, double t){
    if(analytic){
        return function(position, t);
    }
    /*
    Eigen::VectorXd out;
    out.resize(position.size());
    for(int i = 0; i < position.size(); i++ ){

        out(i) = lookuptable[i].Eval(position);
    }
    */
   Eigen::VectorXd F;
   F.resize(2);
   F(0) = lookuptable.Dx(position)*1.41421356237;
   F(1) = lookuptable.Dy(position)*1.41421356237;
   return -F;
}

Eigen::VectorXd Default_Vector(Eigen::VectorXd position, double t){

    return Eigen::VectorXd::Zero(position.size());
    
}
