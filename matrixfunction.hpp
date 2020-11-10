#ifndef MATRIXFUNCTION
#define MATRIXFUNCTION

#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

/*This class takes care of the functions which image is a matrix in BVP*/
class MatrixFunction{

    typedef Eigen::MatrixXd (*pfmatrix)(Eigen::VectorXd , double);

    private:

        //stores the function if it is analytically defined
        pfmatrix function;

        //stores the function if it is stored in a LUT
        std::vector<std::vector<LookUpTable>> lookuptable;

        //True if the function is analytic, false if it is not
        bool analytic;

    public:
    
        //Initializes the class
        MatrixFunction(void);

        /*
        Initialization of the object with the function
        it is suposed to perform.
        input has to be of the kind:
        Eigen::MatrixXd (*input)(Eigen::VectorXd X, double t);
        */
        void Init(pfmatrix input);

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
        Inputs: std::EigenvectorXf, double
        */
        Eigen::MatrixXd Value(Eigen::VectorXd position, 
                    double t);
};  

//Default function which always returns a matrix full of 0.0f
Eigen::MatrixXd Default_Matrix(Eigen::VectorXd position, double t);

#endif 