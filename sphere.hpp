#include <math.h>
#include <eigen3/Eigen/Core>

double Sphere(double* params, 
            Eigen::VectorXd & position, 
            Eigen::VectorXd & exitpoint,
            Eigen::VectorXd & normal);
void Plot_Sphere(double * params, Eigen::VectorXd & position);
