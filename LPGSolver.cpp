#include"LPGSolver.hpp"
LPGSolver::LPGSolver(void){
    h = 0.001;
}
LPGSolver::LPGSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed){
        Configure(boundary_value_problem,boundary_parameters,discretization, seed);
}

void LPGSolver::Update_Stat(double sol_0, double xi)
{

    sol = sol_0 + xi;
    sums[0] += sol;
    sums[1] += sol*sol;
    sums[2] += sol_0;
    sums[3] += sol_0*sol_0;
    sums[4] += xi;
    sums[5] += xi*xi;
    sums[6] += xi*sol_0;
    sums[7] += pow(sol - sol_a,2);
    N_trayectories ++;
    norm = 1.0/N_trayectories;
    err =fabs(mean-sol_a);
    rerr = err/sol_a;
    mean = sums[0] * norm;
    mse = sums[7] * norm;
    var = sums[1]*norm - mean*mean;
    std = sqrt(norm*var);
}

bool LPGSolver::Inside(){
    bool stoppingbc = bvp.boundary.stop(E_P);
    double dist = bvp.boundary.Dist(params, X, E_P, N);

    if(stoppingbc){

      if( dist < -0.0){

          if (t > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
      } else {

        if (t > 0.0) {

              status = stop;

          } else {

              X = E_P;
              status = time_out;

          }

      }

    }else{
      if( dist <= -0.0){

        if (t > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
        } else {
            if (t > 0.0) {
          
              status = reflect;

            } else {

              X = E_P;
              status = time_out;

          }
      }
    }

    switch(status){

      case stop:
        ji_t = 0.0;
        X = E_P;
        return false;
        break;

      case reflect:
        ji_t = dist;
        X = E_P;
        return true;
        break;

      case time_out:
        ji_t = 0.0;
        t = 0.0;
        return false;
        break;

      default:
        ji_t = 0.0;
        return true;
        break;
    }
}
void LPGSolver::LPG_Step(double &rho){
    if (d_k > -rho){
        do{
            Increment_Update();
            N_rngcalls += X.size();
            Xp = X + bvp.b.Value(X,t)*h + sigma*increment;
            omega =  gsl_ran_exponential(rng,2*h); //Exponential distribution with parameter 1/(2*h)
            uc = N.transpose().dot(sigma*increment +bvp.b.Value(X,t)*h);
            nu = 0.5 *(uc+sqrt(pow((N.transpose()*sigma).norm(),2.0)*omega+pow(uc,2.0f)));
            d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
            if (d_k < -0.0) d_k = 0.0;
            ji_t = std::max(0.0,nu+d_k);
            Xp = Xp - ji_t*N;
        }while((Xp - E_P).dot(N)>0.0);
        d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
        if(d_k > 0.0){
            //printf("WARNING: The particle didn't enter in the domain  after Lepingle step.\n");
            Xp = E_Pp;
            ji_t += d_k;
            d_k = 0.0;
        }
    } else {
        do{
            Increment_Update();
            N_rngcalls += X.size();
            Xp = X + bvp.b.Value(X,t)*h + sigma*increment;
            ji_t = 0.0;
            d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
            if(d_k > -0.0) printf("WARNING: Rho value was understimated\n");
        }while(d_k > 0.0);
    }
}

void LPGSolver::RB_Step(double &rho){
    Increment_Update();
    N_rngcalls += X.size();
    Xp = X + bvp.b.Value(X,t)*h + sigma*increment;
    d_p = bvp.boundary.Dist(params, Xp,E_Pp,Np);
    if(d_p >= 0.0)
    {   
        Xp = E_Pp;
        in_t = false;

    }else {
        if (d_k > -rho){
            omega = gsl_rng_uniform(rng);
            nu = exp((-2.0*d_p*d_k)/(h*(N.transpose()*sigma.transpose()).dot(sigma*N)));
            if (omega < nu){
                Xp = E_Pp;
                in_t = false;
            }
        }
        d_k = d_p;
    }
}
double LPGSolver::Solve(Eigen::VectorXd X0, double T_start, double rho, double tolerance){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    N_trayectories = 0;
    N_rngcalls = 0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    Xp.resize(X0.size());
    E_Pp.resize(X0.size());
    Np.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    //Sample test
    for(unsigned int i = 0; i < 50000; i++)
    {   
        Reset(X0, T_start);
        d_k = bvp.boundary.Dist(params,X,E_P,N);
        in_t = true;
        do{ 
            sigma = bvp.sigma.Value(X,t);
            if(bvp.boundary.stop(X)){
                RB_Step(rho);
            }else{
                LPG_Step(rho);
            }
            
            xi += Y*bvp.F.Value(X,t).dot(increment);
            Z += Y*(bvp.f.Value(X,t)*h + bvp.psi.Value(X,N,t)*ji_t);
            Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
            X = Xp;
            N = Np;
            E_P = E_Pp;
            t += - h;
            ji_t = 0.0;
       }while(t > 0.0  && in_t);

       if(in_t){
            sol_0 = Y*bvp.g.Value(X,0.0f)+Z;
       } else {
            sol_0 = Y*bvp.g.Value(X,t)+Z;
       }
       Update_Stat(sol_0,xi);
    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        Reset(X0, T_start);
        d_k = bvp.boundary.Dist(params, X,E_P,N);
        in_t = true;
        do{ 
            sigma = bvp.sigma.Value(X,t);
            if(bvp.boundary.stop(X)){
                RB_Step(rho);
            }else{
                LPG_Step(rho);
            }
            
            xi += Y*bvp.F.Value(X,t).dot(increment);
            Z += Y*(bvp.f.Value(X,t)*h + bvp.psi.Value(X,N,t)*ji_t);
            Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
            X = Xp;
            N = Np;
            E_P = E_Pp;
            t += - h;
            ji_t = 0.0f;
       }while(t > 0.0  && in_t);

       if(in_t){
            sol_0 = Y*bvp.g.Value(X,0.0)+Z;
       } else {
            sol_0 = Y*bvp.g.Value(X,t)+Z;
       }
       Update_Stat(sol_0,xi);
    }while(2.0*std > tolerance*err);
    norm = 1.0/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*norm)*norm;
    var_xi = sums[5]*norm - sums[4]*sums[4]*norm*norm;
    pearson_c = covar/sqrt((sums[3]*norm - (sums[0]-sums[4])*
    (sums[0]-sums[4])*norm*norm)*var_xi);

    return mean;
}

double LPGSolver::Solve(Eigen::VectorXd X0, double rho, double tolerance){
    Solve(X0,INFINITY,rho,tolerance);
    return mean;
}
double LPGSolver::Solve(Eigen::VectorXd X0, double T_start, double rho, unsigned int Ntray){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    Xp.resize(X0.size());
    E_Pp.resize(X0.size());
    Np.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    for(unsigned int i = 0; i < Ntray; i++)
    {
        Reset(X0, T_start);
        d_k = bvp.boundary.Dist(params, X,E_P,N);
        in_t = true;
        do{ 
            sigma = bvp.sigma.Value(X,t);
            if(bvp.boundary.stop(X)){
                RB_Step(rho);
            }else{
                LPG_Step(rho);
            }
            
            xi += Y*bvp.F.Value(X,t).dot(increment);
            Z += Y*(bvp.f.Value(X,t)*h + bvp.psi.Value(X,N,t)*ji_t);
            Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
            X = Xp;
            N = Np;
            E_P = E_Pp;
            t += - h;
            ji_t = 0.0f;
       }while(t > 0.0  && in_t);

       if(in_t){
            sol_0 = Y*bvp.g.Value(X,0.0)+Z;
       } else {
            sol_0 = Y*bvp.g.Value(X,t)+Z;
       }
       Update_Stat(sol_0,xi);
    }
    norm = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*norm)*norm;
    var_xi = sums[5]*norm - sums[4]*sums[4]*norm*norm;
    pearson_c = covar/sqrt((sums[3]*norm - (sums[0]-sums[4])*
    (sums[0]-sums[4])*norm*norm)*var_xi);

    return mean;
}
double LPGSolver::Solve(Eigen::VectorXd X0, double rho, unsigned int Ntray){
  return  Solve(X0, INFINITY, rho, Ntray);
}