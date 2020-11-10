#ifndef BoundaryValueProblem
#define BoundaryValueProblem

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "scalarfunction.hpp"
#include "scalarfunctionN.hpp"
#include "vectorfunction.hpp"
#include "matrixfunction.hpp"
#include "rbfunction.hpp"
#include "boundary.hpp"
typedef double (*pfscalar)(Eigen::VectorXd, double);
typedef double (*pfscalarN)(Eigen::VectorXd, Eigen::VectorXd, double);
typedef Eigen::VectorXd(*pfvector)(Eigen::VectorXd,double);
typedef Eigen::MatrixXd(*pfmatrix)(Eigen::VectorXd, double);
typedef double (*pfbound)(double* params, 
                             Eigen::VectorXd& position, 
                             Eigen::VectorXd& exitpoint,
                             Eigen::VectorXd& normal);
typedef bool (*pfstop)(Eigen::VectorXd);
typedef double (*pRBF)(Eigen::VectorXd, Eigen::VectorXd, double c);
/*
 BVP stores the functions that can be used to compute the solution of 
 such kind of problems using the Feynman-Kak formula.
*/
class BVP
{  

    private:

        std::map<std::string, bool> control;

    public:

        //Functions whose image is scalar s.t. c(X), f(X), u(X), varphi(X)
        ScalarFunction f,c, u, g, p;
        ScalarFunctionN psi, varphi;
        //Functions whose image is a vector s.t b(X), F(x), g(X), psi(X)
        VectorFunction F, mu, b;

        //sigma is the /sigma matrix
        MatrixFunction sigma;
        
        //The boundary is stored in boundary
        Boundary boundary;

        //RBP function for the meshless method
        RBFunction rbf;
        //All the entries of the maps are set as Default. bvplat remains empty
        BVP(void);

        //The object is properly initialized given the functions or the directorys where the look up tables are stored.
        void BVP_init(int dim,
                    std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfscalarN> map_fscalarN,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> map_fmatrix,
                    std::map<std::string, std::string> map_lut);

        void BVP_init(int dim,
                    std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfscalarN> map_fscalarN,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> map_fmatrix,
                    std::map<std::string, std::string> map_lut,
                    pRBF rbfunc);
        
        /*Initialization of the surface variable with analytic boundary
        -boundary is a function s.t. (*double)(double* params, Eigen::VectorXd& position, 
                             Eigen::VectorXd& exitpoint, Eigen::VectorXd& normal)
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXd)
        */
        void Boundary_init(pfbound boundary, pfstop stopf);
        
        /*Initialization of the surface variable with LUT boundary
        -dim is the dimension of the problem
        -boundary is the directory where the lup of the distance is storedAh per
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXd)
        */
        void Boundary_init(int dim, std::string boundary , pfstop stopf);

};




#endif