#ifndef VECTORFUNCTION
#define VECTORFUNCTION

#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the functions which image is a vector in BVP*/
class VectorFunction{

    typedef Eigen::VectorXd (*pfvector)(Eigen::VectorXd,double);

    private:

        //stores the function if it is analytically defined
        pfvector function;

        //stores the function if it is stored in a LUT
        /*std::vector<LookUpTable> lookuptable;*/
        LookUpTable lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:

        //Initializes the class
        VectorFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        Eigen::VectorXd (*input)(Eigen::VectorXd X, Eigen::VectorXd N);
        */
        void Init(pfvector input);

        /*
        Initialization of the object with the lut where
        The values of the function are stored
        dim has to be an unsigned integer
        input has to be a std::string
        */
        void Init(int dim, 
                  std::string input);

        /*
        Returns the value of the function in X with
        normal vector N.
        Inputs: std::EigenvectorXf, std::EigenvectorXd
        */
        Eigen::VectorXd Value(Eigen::VectorXd position,double t);
};  

//Default function which always returns a vector full of 0.0f
Eigen::VectorXd Default_Vector(Eigen::VectorXd position, double t);

#endif 