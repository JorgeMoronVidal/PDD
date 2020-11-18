#include "FKACSolver.hpp"
#include "stencil.hpp"
class GMSolver: public EMFKAC {
    private:
        /*-dist: Distance to the boundary*/
        double dist, norm;
        /*-Are bc's stopping?*/
        bool stoppingbc;
        bool Inside();
        /*Updates statistical quantities*/
    void Update_Stat(double sol_0, double xi);
    public:
        GMSolver(void);
        /*Class initialization
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process 
        -seed is the RNG seed*/  
        GMSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator till 2*std < tolerance*err*/
        double Solve(Eigen::VectorXd X0, double T_start, double tolerance);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, double tolerance);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator with N trayectories*/
        double Solve(Eigen::VectorXd X0, double T_start, unsigned int Ntray);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, unsigned int Ntray);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        void Solve(Eigen::VectorXd X0, double c2, Stencil stencil, 
                    std::vector<int> & G_j, std::vector<double> & G,
                    double & B, unsigned int Ntray);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator with N trayectories*/
        void Solve(Eigen::VectorXd X0, double T_start, double c2,
                    Stencil & stencil, std::vector<int> & G_j, 
                    std::vector<double> & G, double &B,  unsigned int Ntray);
        /*Convergence test for a given parabolic BVP*/
        void Test(std::string filename, Eigen::VectorXd X0, double T_start, 
                  double tolerance, double h0, unsigned int Nsamples);
        /*Convergence test for a given elliptic BVP*/
        void Test(std::string filename, Eigen::VectorXd X0, double tolerance, 
                  double h0, unsigned int Nsamples);
};