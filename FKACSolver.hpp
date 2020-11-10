#ifndef EulerMayuramaFeynmanKac
#define EulerMaruyamaFeynmanKac
#include <stdlib.h>
#include <iostream> 
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "BVP.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__);}
enum outcome {in, stop, reflect, time_out};
class EMFKAC
{
protected:
    /*Boundary value problem's equations'*/
    BVP bvp;
    /*-params: Parameters of the surface
      -t:Time (For parabolic equations).
      -h: Time discretization.
      -sqrth: square root of h.
      -Y: FKAK Y.
      -Z: FKAK Z.
      -ji_t auxiliary time in order to deal with reflecint BC's.
      -xi: Control variates variable.
      -params stores the parameters of the surface
      -sums stores important quantities to compute statistical measurements
      -sola_a stores the analytic solution if available.*/
    double t, h, sqrth, Y, Z, ji_t, xi, *params, sums[10], sol, sol_0;
    /*-N: Normal vector to the boundary.
      -X: Position.
      -Variance reduction mu function.
      -increment: Random increment in position.*/
    Eigen::VectorXd N, E_P, X, mu, increment;
    /*Sigma matrix*/
    Eigen::MatrixXd sigma;
    /*RNG type*/
    const gsl_rng_type * T;
    /*RNG*/
    gsl_rng * rng;
    /*Controls the status of the particle*/
    outcome status;
    /*Updates increment with random numbers*/
    void Increment_Update(void);
    /*Reset of FKAK */
    void Reset(Eigen::VectorXd X0, double T_start);
    /*One plain step of the Euler's algorithm*/
    void Step(void);
    /*One step of the Euler's discretization with Variance Reduction*/
    void VR_Step(void);
    /*One step of the Euler's discretization with Control Variates*/
    void CV_Step(void);
    /*One step of the Euler's discretization with Variance Reduction and 
    Control Variates*/
    void VR_CV_Step(void);
public:
    /*-mean and solution of the solver
      -variance of the solution 
      -standard deviation of the solution 
      -absolute error (with sign)
      -relative error (with sign)
      -mean square error */
    double mean, var, std, err, rerr, mse, covar, var_xi, pearson_c, sol_a;
    /*-Number of calls to the random number generator.
      -Number of trayectories computed.
    */
    unsigned int N_rngcalls, N_trayectories;
    /*Empty initizialization*/
    EMFKAC(void);
    /*Proper initialization with default random seed equal to 0
    -boundary_value_problem is a BVP object which stores all problem's equations
    -surface_parameters stores the parameters for the boundary construction
    -discretization stores the time discretization for the Stochastic process 
    -seed is the random number generator seed*/
    EMFKAC(BVP boundary_value_problem,
      std::vector<double> surface_parameters,
      double discretization,
      unsigned int seed);
    /*Changes class problem and parameters 
    -boundary_value_problem is a BVP object which stores all problem's equations
    -boundary_parameters stores the parameters for the boundary construction
    -discretization stores the time discretization for the Stochastic process
    -seed is the random number generator seed*/
    void Configure(BVP boundary_value_problem,
      std::vector<double> surface_parameters,
      double discretization,
      unsigned int seed);
    //~EMFKAC(); 
};
#endif