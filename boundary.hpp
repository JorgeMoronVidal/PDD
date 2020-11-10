#ifndef BOUND
#define BOUND

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
/*Boundary is a class which stores the boundary-related functions*/
class Boundary{

    typedef double (*pfbound)(double*, Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
    typedef bool (*pfstop)(Eigen::VectorXd);

    private:

        /*If the distance function is given by an analytic 
        function, it is stored in Distance_Analytic */
        pfbound Distance_Analytic;

        /*If the distance function is given by a look up 
        table, it is stored in Distance_Numeric */
        LookUpTable Distance_Numeric;
        
        /*True if Distance_Analytic gives the distance
          False if Distance_Numeric gives the distance*/
        bool analytic;

    public:
        //True if stopping, false if reflecting
        pfstop stop;

        //Initialization by default
        Boundary(void);

        //Initialization with Analytic Distance function
        void _init_(pfbound fbound, pfstop fstop);

        //Initialization Distance stored in a look up table
        void _init_(int dim, std::string lupbound, pfstop fstop);

        /*-Returns distance to the boundary (Negative if inside, positive if outside)
          -Modifies the normal vector
          -Modifies exitpoint
          -If outside, particle is placed on the boundary again*/
        double Dist(double* params, 
                   Eigen::VectorXd & position, 
                   Eigen::VectorXd & exitpoint,
                   Eigen::VectorXd & normal);

};
#endif