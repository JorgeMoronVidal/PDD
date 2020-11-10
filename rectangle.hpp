#ifndef RECTANGLE
#define RECTANGLE
#include <math.h>
#include <eigen3/Eigen/Core>
#include <vector>
#include <iostream>
double Rectangle2D(double* params, 
            Eigen::VectorXd & position, 
            Eigen::VectorXd & exitpoint,
            Eigen::VectorXd & normal);
bool Stopping(Eigen::VectorXd position);

#endif