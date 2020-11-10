#ifndef BORDERHEADERDEF
#define BORDERHEADERDEF
#include <eigen3/Eigen/Core>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/*2D look up table class which uses gsl spline2D*/
class LookUpTable{

    private:
    //Stores  the grid points
    double** x;
    //Stores the value of the function on the grid
    double* z;
    //Lengt of the grid in each direction
    int *len;
    //gsl spline instances to interpolate the grid values
    const gsl_interp2d_type *T;
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    
    public:
    
    //Initialization of the clase
    LookUpTable(void);
    /*-dim is the dimension of the problem [2, at this point in the developement]
      -file is the file where the look up tables are stored*/
    void Init(int dim, std::string directory);
    /*
    Given a position or a position, an exitpoint and a normal 
    vector, it 
    */
    double Dx(Eigen::VectorXd position);
    double Dy(Eigen::VectorXd position);
    double Eval(Eigen::VectorXd position);
    double Eval(double* params, 
               Eigen::VectorXd & position, 
               Eigen::VectorXd & exitpoint,
               Eigen::VectorXd & normal);
};
#endif