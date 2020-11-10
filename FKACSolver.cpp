#include "FKACSolver.hpp"
//#define PLOT
EMFKAC::EMFKAC()
{
    //Default values
    h = 0.001;
    sqrth = sqrt(h);
    t = 0.0;
    T = gsl_rng_default; //Mt1997
    rng = gsl_rng_alloc(T);
    N_rngcalls = 0;
    N_trayectories = 0;
    T = gsl_rng_default; //Mt1997
    rng = gsl_rng_alloc(T);
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
}

EMFKAC::EMFKAC(BVP boundary_value_problem,
      std::vector<double> surface_parameters,
      double discretization,
      unsigned int seed)
{
    Configure(boundary_value_problem, surface_parameters, discretization, seed);
}


void EMFKAC::Configure(BVP boundary_value_problem,
      std::vector<double> surface_parameters,
      double discretization,
      unsigned int seed)
{  
    h = discretization;
    sqrth = sqrt(h);
    bvp = boundary_value_problem;
    //delete params;
    params = new double[surface_parameters.size()];
    for(unsigned int i = 0; i < surface_parameters.size(); i++) params[i] = surface_parameters[i];
    t = INFINITY;
    N_rngcalls = 0;
    N_trayectories = 0;
    T = gsl_rng_default; //Mt1997
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed);
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
}

void EMFKAC::Increment_Update(void){
    for(int j = 0; j < increment.size(); j++)
    {
        increment(j) = gsl_ran_gaussian_ziggurat(rng,1)*sqrth;
    }
}


void EMFKAC::Reset(Eigen::VectorXd X0, double T_start)
{
    X = X0;
    t = T_start;
    Y = 1.0;
    Z = 0.0;
    xi = 0.0;
    ji_t = 0.0;
}

void EMFKAC::Step(void){
    
    sigma = bvp.sigma.Value(X,t);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += bvp.b.Value(X,t)*h + sigma*increment;
    t -= h;
}

void EMFKAC::VR_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    mu = bvp.mu.Value(X,t);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + Y*bvp.mu.Value(X,t).transpose()*increment + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += (bvp.b.Value(X,t) - sigma*mu)*h + sigma*increment;
    t -= h;
}

void EMFKAC::CV_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    xi += Y*bvp.F.Value(X,t).dot(increment);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)*Y*ji_t;
    Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += bvp.b.Value(X,t)*h + sigma*increment;
    t -= h;
}

void EMFKAC::VR_CV_Step(void)
{   
    sigma = bvp.sigma.Value(X,t);
    mu = bvp.mu.Value(X,t);
    xi += Y*bvp.F.Value(X,t).dot(increment);
    Z += bvp.f.Value(X,t)*Y*h + bvp.psi.Value(X,N,t)* Y *ji_t;
    Y += bvp.c.Value(X,t)*Y*h + Y*mu.transpose()*increment + bvp.varphi.Value(X,N,t) * Y * ji_t;
    X += (bvp.b.Value(X,t) - sigma*mu)*h + sigma*increment;
    t -= h;
}

