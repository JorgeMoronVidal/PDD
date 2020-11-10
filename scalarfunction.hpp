#ifndef SCALARFUNCTION
#define SCALARFUNCTION

#include <iostream>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the scalar functions in BVP*/
class ScalarFunction{

    typedef double (*pfscalar)(Eigen::VectorXd , double );

    private:

        //stores the function if it is analytically defined
        pfscalar function;

        //stores the function if it is stored in a LUT
        LookUpTable lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:

        //Initializes the class
        ScalarFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        double (*input)(Eigen::VectorXd X, double);
        */
        void Init(pfscalar input);

        /*
        Initialization of the object with the lut where
        The values of the function are stored
        dim has to be an integer
        input has to be a std::string
        */
        void Init(int dim, 
                  std::string input);

        /*
        Returns the value of the function in X
        Inputs: std::EigenvectorXd, double 
        */
        double Value(Eigen::VectorXd position, 
                    double t);
};  
//Default function which always returns  0.0f
double Default_Scalar(Eigen::VectorXd X,double t);
#endif 

