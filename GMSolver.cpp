#include"GMSolver.hpp"
#include"rectangle.hpp"
GMSolver::GMSolver(void){
    h = 0.001;
}
GMSolver::GMSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed){
        Configure(boundary_value_problem,boundary_parameters,discretization, seed);
}

void GMSolver::Update_Stat(double sol_0, double xi)
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

bool GMSolver::Inside(){
    bool stoppingbc = bvp.boundary.stop(E_P);
    double dist = bvp.boundary.Dist(params, X, E_P, N);

    if(stoppingbc){

      if( dist < -0.5826*(N.transpose()*sigma).norm()*sqrth){

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
      if( dist <= -0.5826*(N.transpose()*sigma).norm()*sqrth){

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
double GMSolver::Solve(Eigen::VectorXd X0, double T_start, double tolerance){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    N_trayectories = 0;
    N_rngcalls = 0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        Reset(X0, T_start);
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        Reset(X0, T_start);
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());
        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
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

double GMSolver::Solve(Eigen::VectorXd X0, double tolerance){
    Solve(X0,INFINITY,tolerance);
    return mean;
}
double GMSolver::Solve(Eigen::VectorXd X0, double T_start, unsigned int Ntray){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    for(unsigned int i = 0; i < Ntray; i++)
    {
        Reset(X0, T_start);

        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
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
double GMSolver::Solve(Eigen::VectorXd X0, unsigned int Ntray){
  return  Solve(X0, INFINITY, Ntray);
}

void GMSolver::Solve(Eigen::VectorXd X0, double c2, Stencil stencil, 
             std::vector<int> & G_j, std::vector<double> & G, double &B, 
             unsigned int Ntray){
Solve(X0, INFINITY, c2, stencil, G_j, G, B, Ntray);
}
void GMSolver::Solve(Eigen::VectorXd X0, double T_start, double c2,
             Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G, double & B,  unsigned int Ntray){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    //Boundary of the stencil
    Boundary sten_boundary;
    sten_boundary._init_(Rectangle2D, Stopping);
    //G and B are emptied
    G_j.clear();
    G.clear();
    double b;
    B = 0.0;
    stencil.Reset();
    /*Control variable 
    -0 if trayectory still inside
    -1 if trayectory exits by stencil's boundary
    -2 if trayectory exits by problem's boundary*/
    int16_t control;
    Reset(X0,T_start);
    if(bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N) >= -1E-06){
                B = bvp.g.Value(X0,t);
    } else{
      for(int i = 0; i < Ntray; i++){
        Reset(X0,T_start);
        control = 0;
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
            if(t <= 0){
                control = 3;
            }
            if(bvp.boundary.Dist(params, X, E_P, N) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                control = 2;
                break;
            }
            if(sten_boundary.Dist(stencil.stencil_parameters, X, E_P, N) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                control = 1;
                break;
            }
        } while(control == 0);
        
        switch (control) {        
            case 1: 
                stencil.G_update(E_P,Y,bvp,c2);
                b += Z + xi;
            break;        
            case 2:
                b += Z + xi + Y*bvp.g.Value(E_P,t);
            break;
            case 3:
                b += Z + xi + Y*bvp.p.Value(X,t);
            break;
            default : 
                std::cout << "Something went wrong while solving";
        }
      }
      stencil.G_return_withrep(G_j, G, Ntray);
      B = b/Ntray;
  }
}


void GMSolver::Test(std::string filename, Eigen::VectorXd X0, double T_start, 
                  double tolerance, double h0,  unsigned int Nsamples){
    FILE *fp;
    fp = fopen(filename.c_str(),"w");
    time_t now = time(NULL);
    char *date = ctime(&now);
    int len = strlen(date);
    date[len-1] = ' ';
    printf("\n");
    printf("**********************************************************\n");
    printf("*** FKAKsolver class file version 0.1    produced by J. Moron ***\n");
    printf("*** C++ elliptic_test on %s         ***\n",date);
    printf("**********************************************************\n");
    printf("\n");
    printf("**********************************************************\n");
    printf("*** Convergence tests ************************************\n");
    printf("****** ****************************************************\n");
    printf("\n    h \t   exa_sol  num_sol   variance   pcoef \t    std       error   ");
    printf(" r.error    mse       cost \n---------------------------------");
    printf("-------------------------------------------------------------\n");
    fprintf(fp,"h,exact_sol,num_sol,variance,pcoef,std,err,rerr,mse,cost\n");
    fclose(fp);
    for(unsigned int i = 0; i < Nsamples; i++){
        h = h0 * pow(2,-1.0*i);
        sqrth = sqrt(h);
        Solve(X0, T_start, tolerance);
        fopen(filename.c_str(),"a");
        printf("%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e \n",
        h,sol_a, mean,var,pearson_c,std,err,rerr,mse,1.0*N_rngcalls);
        fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
        h,sol_a, mean,var,pearson_c,std,err,rerr,mse,1.0*N_rngcalls);
        fclose(fp);
    }
}

void GMSolver::Test(std::string filename, Eigen::VectorXd X0, double tolerance, 
                  double h0, unsigned int Nsamples){
        Test(filename,X0,INFINITY,tolerance,h0,Nsamples);
    }




