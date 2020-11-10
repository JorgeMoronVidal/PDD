#include "FKACSolver.hpp"
#define NOISYBOUNDARY
//#define NOISYPSEUDOSPLIT
#include <ctime>
class MLVLSolver: public EMFKAC {
    private:
    /*-dist: Distance to the boundary*/
    double dist, norm;
    /*-Are bc's stopping?*/
    bool stoppingbc;
    bool Inside(Eigen::VectorXd &position,
        Eigen::VectorXd &normal, 
        Eigen::VectorXd &exitpoint,
        double &time,
        double &ji_time, 
        double sqrtdiscretization);
    /*Multilevel related variables and functions*/
    /*-l level
        -option 1 -- sub-sampling, no offset
        -2 -- nothing
        -3 -- sub-sampling, offset
      -M is the number of fine increments that sums up acoarse increment.
    */
    uint16_t l, option, M, D, L;
    uint32_t dNl[21], Nl[21];
    Eigen::VectorXd Xc, incrementc, Nc, E_Pc;
    double h0, hc, sqrthc, Yc, Zc, xic, tc, ji_tc, Pf, Pc, dP, Cl[21],ml[21], Vl[21], NlCl[21], mse[21];
    /*Updates coarsed increment with random numbers*/
    void Incrementc_Update(void);
    /*Step of FKAK formula for corased set of variables*/
    void Stepc(void);
    /*Reset of FKAK */
    void Resetc(Eigen::VectorXd X0);
    void Resetc(Eigen::VectorXd X0, double t_start);
    /*Computes one multilevel trayectory of level l for elliptic equations*/
    void Solve_l(Eigen::VectorXd X0);
    /*Computes one multilevel trayectory of level l for parabolic equations*/
    void Solve_l(Eigen::VectorXd X0, double t_start);
    void Solve_l(Eigen::VectorXd X0, double t_start, double noise_mean, double noise_var);
    /*Linear regression routine*/
    void Regression(unsigned int N_sample, double *x, double *y, double &a, double &b);

    public:
        MLVLSolver(void);
        /*Class initialization
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process
        -M_factor is the Giles multilevel M 
        -seed is the RNG seed*/  
        MLVLSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed);
        MLVLSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double initial_discretization,
            uint16_t M_factor,
            unsigned int seed);
        /*Class configuration
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process
        -M_factor is the M variable from Giles 2008 paper [optimal theoretical value is 4].
        -seed is the RNG seed
        */
        void Init(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double initial_discretization,
            uint16_t M_factor,
            unsigned int seed);
        /*Solves an elliptic equation on inital point X0 using Multilevel 
        algorithm*/
        double Solve(Eigen::VectorXd X0, int Lmin, int Lmax, int N0, double eps,
                    double alpha_0,double beta_0,double gamma_0);
        /*Solves a parabollic equation on inital point X0 using Multilevel 
        algorithm*/
        double Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
                    double alpha_0,double beta_0,double gamma_0);
        double Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
                    double alpha_0,double beta_0,double gamma_0,std::string fn_results, std::string fn_nl);
        double Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
                    double noise_mean, double noise_variance, double alpha_0,double beta_0,double gamma_0,
                    std::string fn_results, std::string fn_nl);
        /*Multilevel test for the elliptic solver*/
        void Test(Eigen::VectorXd X0, uint32_t N_test, uint32_t L_test, uint32_t N0, double *Eps, uint32_t Lmin, uint32_t Lmax, std::string fname);
        /*Multilevel test for the parabolic solver*/
        void Test(Eigen::VectorXd X0, double t_start, uint32_t N_test,uint32_t L_test, uint32_t N0, double *Eps, uint32_t Lmin, uint32_t Lmax, std::string fname);
};