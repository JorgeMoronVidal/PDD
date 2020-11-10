#include "FKACSolver.hpp"

class LPGSolver: public EMFKAC {
    private:
        /*-dist: Distance to the boundary*/
        double uc,nu,omega,d_p,d_k, dist, norm;
        /*-Are bc's stopping?*/
        bool stoppingbc, in_t;
        Eigen::VectorXd Xp,E_Pp,Np;
        bool Inside();
        /*Updates statistical quantities*/
        void Update_Stat(double sol_0, double xi);
        /*Step of Lepingle stimator*/
        void LPG_Step(double &rho);
        /*Step of the Random Bridge stimator*/
        void RB_Step(double &rho);

    public:
        LPGSolver(void);
        /*Class initialization
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process 
        -seed is the RNG seed*/  
        LPGSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator till 2*std < tolerance*err*/
        double Solve(Eigen::VectorXd X0, double T_start, double rho, double tolerance);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, double rho, double tolerance);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator with N trayectories*/
        double Solve(Eigen::VectorXd X0,  double T_start, double rho,  unsigned int Ntray);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, double rho, unsigned int Ntray);
};