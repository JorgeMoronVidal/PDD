#include"MLVLSolver.hpp"


//#define PLOT
MLVLSolver::MLVLSolver(void){
    h = 0.001;
}
MLVLSolver::MLVLSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double initial_discretization,
            unsigned int seed){
    M = 4;
    h0 = initial_discretization;
    Configure(boundary_value_problem,boundary_parameters,initial_discretization, seed);
}
MLVLSolver::MLVLSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double initial_discretization,
            uint16_t M_factor,
            unsigned int seed){
    Init(boundary_value_problem, boundary_parameters, initial_discretization, M_factor,seed);
}
void MLVLSolver::Init(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double initial_discretization,
            uint16_t M_factor,
            unsigned int seed){
    h0 = initial_discretization;
    M = M_factor;
    Configure(boundary_value_problem,boundary_parameters,initial_discretization, seed);
}

void MLVLSolver::Stepc(void){

    sigma = bvp.sigma.Value(Xc,tc);
    mu = bvp.mu.Value(Xc,tc);
    xic += Yc*bvp.F.Value(Xc,tc).dot(incrementc);
    Zc += bvp.f.Value(Xc,tc)*Yc*hc + bvp.psi.Value(Xc,Nc,tc)*Yc*ji_tc;
    Yc += bvp.c.Value(Xc,tc)*Yc*hc + Yc*mu.transpose()*incrementc + bvp.varphi.Value(Xc,Nc,t) * Yc * ji_tc;
    Xc += (bvp.b.Value(Xc,tc) - sigma*mu)*hc + sigma*incrementc;
    tc -= hc;

}

bool MLVLSolver::Inside(Eigen::VectorXd &position,
              Eigen::VectorXd &normal,
              Eigen::VectorXd &exitpoint,
              double &time,
              double &ji_time, 
              double sqrtdiscretization){
    stoppingbc = bvp.boundary.stop(exitpoint);
    dist = bvp.boundary.Dist(params, position, exitpoint, normal);
    
    if(stoppingbc){

    if( dist < -0.5826*(normal.transpose()*sigma).norm()*sqrtdiscretization){
      //if( dist < -0.0){

          if (time > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
      } else {

        if (time > 0.0) {

              status = stop;

          } else {

              position = exitpoint;
              status = time_out;

          }

      }

    }else{
      if( dist <= -0.5826*(normal.transpose()*sigma).norm()*sqrtdiscretization){

        if (time > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
        } else {
            if (time > 0.0) {
          
              status = reflect;

            } else {

              position = exitpoint;
              status = time_out;

          }
      }
    }

    switch(status){

      case stop:
        ji_time = 0.0;
        position = exitpoint;
        return false;
        break;

      case reflect:
        ji_time = dist;
        position = exitpoint;
        return true;
        break;

      case time_out:
        ji_time = 0.0;
        time = 0.0;
        return false;
        break;

      default:
        ji_time = 0.0;
        return true;
        break;
    }
}

void MLVLSolver::Incrementc_Update(void){
    for(int j = 0; j < incrementc.size(); j++)
    {
        incrementc(j) = gsl_ran_gaussian_ziggurat(rng,1)*sqrthc;
    }
}


void MLVLSolver::Resetc(Eigen::VectorXd X0){
  Resetc(X0,INFINITY);
}

void MLVLSolver::Resetc(Eigen::VectorXd X0, double T_start){
  xic = 0.0;
  Xc = X0;
  Yc = 1.0;
  Zc = 0.0;
  ji_tc = 0.0;
  tc = T_start;
}

void MLVLSolver::Solve_l(Eigen::VectorXd X0){
  Solve_l(X0,INFINITY);
}

void MLVLSolver::Solve_l(Eigen::VectorXd X0, double t_start){
  int16_t split, in_f, in_c, set_f, set_c;
  uint32_t counter = 0;
  double t_split, Y_split, Z_split, xi_split, split_c, split_f, N_steps;
  Eigen::VectorXd X_split;
  increment.resize(X0.size());
  incrementc.resize(X0.size());
  N.resize(X0.size());
  E_P.resize(X0.size());
  Nc.resize(X0.size());
  E_Pc.resize(X0.size());
  D = X0.size();
  h = h0 / pow(M,l);
  sqrth = sqrt(h);
  hc = h0 / pow(M,(l-1));
  sqrthc = sqrt(hc);
  //printf(" **** l, Nl = %u, %u **** h = %f hc = %f\n",l,dNl[l],h,hc);
  sol_a = bvp.u.Value(X0,t_start);
  sol_0 = 0.0;
  if (l == 0){
      sol_0 = 0;
  } else {
      for(int itl = 0; itl < l-1; itl++) {
        sol_0 += ml[itl];
      }
  }
  option = 1;
  for (unsigned int k=0; k<10; k++) sums[k] = 0.0;
  if (option==1 || option==3) {
    //split = 4 * (1<<l);
    split = 1<<l;
    //split = 1;
  }
  else if (option==2) {
    split = 1;
  }
  for (unsigned int np = 0; np < dNl[l]; np++) {
    in_f = 1;
    in_c = 1;
    set_f = 0;
    set_c = 0;
    Pf = 0.0;
    Pc = 0.0;
    counter = 0;
    split_c = 0.0;
    split_f = 0.0;
    N_steps = 0;
    Reset(X0, t_start);
    Resetc(X0, t_start);
    // level 0
    if (l==0) {
      do {
        Increment_Update();
        counter += D;
        VR_CV_Step();
        N_steps ++;
        if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
      } while (in_f);
      if(status == time_out){
          Pf = Y*bvp.p.Value(X,0.0)+Z+xi;
      } else {
          Pf = Y*bvp.g.Value(E_P,t)+Z+xi;
      }
    }
    // level l>0
    else { 
      do {
        for (int d=0; d<D; d++) incrementc(d) = 0.0;
        for (int m=0; m<M; m++){
          counter += D;
          Increment_Update();
          incrementc += increment;
          VR_CV_Step();
          N_steps ++;
          if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
          if ( (!in_f) && (!set_f) ) {
            set_f = 1;
            if(status == time_out){
                Pf = Y*bvp.p.Value(X,0.0)+Z;
            } else {
                Pf = Y*bvp.g.Value(E_P,t)+Z;
            }
            break;
          }
      	}
        Stepc();
        if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
        if ( (!in_c) && (!set_c) ) {
          set_c = 1;
          if(status == time_out){
                Pc = Yc*bvp.p.Value(Xc,0.0)+Zc;
          } else {
                Pc = Yc*bvp.g.Value(E_Pc,tc)+Zc;
          }
        }
      } while (in_f && in_c);
      // split continuation paths
      if (in_f) {
        t_split = t;
        X_split = X;
        Y_split = Y;
        Z_split = Z;
        xi_split = xi;
        split_f = (double) split;
        for (int s=0; s<split; s++) {
	        // reset state at split
          in_f = 1;
          Reset(X_split, t_split);
          Y = Y_split;
          Z = Z_split;
          xi = xi_split;
	        // continue until exit
          do {
            counter += D;
            Increment_Update();
            VR_CV_Step();
            N_steps+=1/split;
            if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
          } while (in_f);
          if(status == time_out){
                Pf += (Y*bvp.p.Value(X,0.0)+Z+xi)/((double) split);
          } else {
                Pf += (Y*bvp.g.Value(E_P,t)+Z+xi)/((double) split);
          }

        }
      }

      if (in_c) {
        t_split = tc;
        X_split = Xc;
        Y_split = Yc;
        Z_split = Zc;
        xi_split = xic;
        split_c = (double) split;
        for (int s=0; s<split; s++) {
	        // reset state at split
          in_c = 1;
          Resetc(X_split, t_split);
          Yc = Y_split;
          Zc = Z_split;
          xic = xi_split;
	        // continue until exit
          do {
            counter += D;
            Incrementc_Update();

            Stepc();
            if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
          } while (in_c);
          if(status == time_out){
                Pc += (Yc*bvp.p.Value(Xc,0.0)+Zc+xic)/ ((double) split);
          } else {
                Pc += (Yc*bvp.g.Value(E_Pc,tc)+Zc+xic)/((double) split);

          }
        }
      }
    }
    //printf("Split c %f Split f %f \n", split_c, split_f);
    dP = Pf-Pc;
    sol = sol_0 + dP;
    sums[0] += counter; // add number of RNG calls
    sums[1] += dP;
    sums[2] += dP*dP;
    sums[3] += dP*dP*dP;
    sums[4] += dP*dP*dP*dP;
    sums[5] += split_f;
    sums[6] += split_c;
    sums[7] += pow(sol-sol_a,2.0); //MSE
    sums[8] += N_steps;
  }
}
/*void MLVLSolver::Solve_l(Eigen::VectorXd X0, double t_start){
  int16_t split, in_f, in_c, set_f, set_c;
  uint32_t counter = 0;
  double t_split, Y_split, Z_split, xi_split, split_c, split_f, N_steps, summ_vr,
  sqsumm_vr, sqsumm_nvr, summ_nvr, var_vr, var_nvr;
  Eigen::VectorXd X_split;
  increment.resize(X0.size());
  incrementc.resize(X0.size());
  N.resize(X0.size());
  E_P.resize(X0.size());
  Nc.resize(X0.size());
  E_Pc.resize(X0.size());
  D = X0.size();
  h = h0 / pow(M,l);
  sqrth = sqrt(h);
  hc = h0 / pow(M,(l-1));
  sqrthc = sqrt(hc);
  //printf(" **** l, Nl = %u, %u **** h = %f hc = %f\n",l,dNl[l],h,hc);
  sol_a = bvp.u.Value(X0,t_start);
  sol_0 = 0.0;
  if (l == 0){
      sol_0 = 0;
  } else {
      for(int itl = 0; itl < l-1; itl++) {
        sol_0 += ml[itl];
      }
  }
  option = 1;
  for (unsigned int k=0; k<10; k++) sums[k] = 0.0;
  if (option==1 || option==3) {
    //split = 4 * (1<<l);
    split = 1<<l;
    //split = 1;
  }
  else if (option==2) {
    split = 1;
  }
  #ifdef PLOT
  char filename[256];
  FILE *file_plot;
  unsigned int j;
  #endif
  for (unsigned int np = 0; np < dNl[l]; np++) {
    #ifdef PLOT
    if(np < 10){
      sprintf(filename,"Output/trayectories/trayectory_MLVL.csv");
      file_plot = fopen(filename, "w");
      for(j = 0; j < X0.size(); j++) fprintf(file_plot,"Xf_%u,E_Pf_%u,Nf_%u,",j,j,j);
      for(j = 0; j < X0.size(); j++) fprintf(file_plot,"Xc_%u,E_Pc_%u,Nc_%u,",j,j,j);
      fprintf(file_plot,"\n");
      fclose(file_plot);
    }
    #endif 
    in_f = 1;
    in_c = 1;
    set_f = 0;
    set_c = 0;
    Pf = 0.0;
    Pc = 0.0;
    counter = 0;
    split_c = 0.0;
    split_f = 0.0;
    N_steps = 0;
    // level 0
    Reset(X0, t_start);
    Resetc(X0, t_start);
    #ifdef PLOT
    if(np < 10){
      file_plot = fopen(filename, "a");
      for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
      for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
      fprintf(file_plot,"\n");
      fclose(file_plot);
    }
    #endif 
    if (l==0) {
      do {
        Increment_Update();
        counter += D;
        VR_CV_Step();
        N_steps ++;
        #ifdef PLOT
        if(np < 10){
          file_plot = fopen(filename, "a");
          for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
          for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
          fprintf(file_plot,"\n");
          fclose(file_plot);
        }
        #endif 
        if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
      } while (in_f);
      if(status == time_out){
          Pf = Y*bvp.p.Value(X,0.0)+Z+xi;
      } else {
          Pf = Y*bvp.g.Value(E_P,t)+Z+xi;
      }
      #ifdef PLOT
      if(np < 10){
        file_plot = fopen(filename, "a");
        for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
        for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
        fprintf(file_plot,"\n");
        fclose(file_plot);
      }
      #endif
    }

      // level l>0

    else {

      Reset(X0, t_start);
      Resetc(X0, t_start);
      #ifdef PLOT
      if(np < 10){
        file_plot = fopen(filename, "a");
        for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
        for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
        fprintf(file_plot,"\n");
        fclose(file_plot);
      }
      #endif 
      do {

        for (int d=0; d<D; d++) incrementc(d) = 0.0;
        for (int m=0; m<M; m++){
          counter += D;
          Increment_Update();
          incrementc += increment;
          VR_CV_Step();
          N_steps ++;
          #ifdef PLOT
          if(np < 10){
            file_plot = fopen(filename, "a");
            for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
            for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
            fprintf(file_plot,"\n");
            fclose(file_plot);
          }
          #endif 
          if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
          if ( (!in_f) && (!set_f) ) {
            set_f = 1;
            if(status == time_out){
                Pf = Y*bvp.p.Value(X,0.0)+Z;
            } else {
                Pf = Y*bvp.g.Value(E_P,t)+Z;
            }
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif 
            break;
          }
      	}
        Stepc();
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
        if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
        if ( (!in_c) && (!set_c) ) {
          set_c = 1;
          if(status == time_out){
                Pc = Yc*bvp.p.Value(Xc,0.0)+Zc;
          } else {
                Pc = Yc*bvp.g.Value(E_Pc,tc)+Zc;
          }
          #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,",,,");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
          #endif
        }
      } while (in_f && in_c);
      // split continuation paths 
      if (in_f) {
        t_split = t;
        X_split = X;
        Y_split = Y;
        Z_split = Z;
        xi_split = xi;
        if(l < 4){
          split_f = (double) split;
          for (int s=0; s<split; s++) {
            #ifdef PLOT
            if(np < 10){
              sprintf(filename,"Output/trayectories/split_%d.csv",s);
              file_plot = fopen(filename, "w");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"X_%u,E_P_%u,N_%u,",j,j,j);
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // reset state at split
            in_f = 1;
            Reset(X_split, t_split);
            Y = Y_split;
            Z = Z_split;
            xi = 0.0;
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // continue until exit
            do {
              counter += D;
              Increment_Update();
              VR_CV_Step();
              N_steps+=1/split;
              #ifdef PLOT
              if(np < 10){
                file_plot = fopen(filename, "a");
                for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
                fprintf(file_plot,"\n");
                fclose(file_plot);
              }
              #endif
              if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
            } while (in_f);
            if(status == time_out){
              Pf += (Y*bvp.p.Value(X,0.0)+Z+xi)/((double) split);
              sqsumm_vr += pow(Y*bvp.p.Value(X,0.0)+Z+xi,2)/((double) split);
            } else {
              Pf += (Y*bvp.g.Value(E_P,t)+Z+xi)/((double) split);
              sqsumm_vr += pow(Y*bvp.p.Value(E_P,t)+Z+xi,2)/
              ((double) split);
            }
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
          }
          var_vr = Pf-sqsumm_vr;
        } else {
          sqsumm_vr = 0.0;
          sqsumm_nvr = 0.0;
          summ_vr = 0.0;
          summ_nvr = 0.0;
          for (int s=0; s<10; s++) {
            #ifdef PLOT
            if(np < 10){
              sprintf(filename,"Output/trayectories/split_%d.csv",s);
              file_plot = fopen(filename, "w");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"X_%u,E_P_%u,N_%u,",j,j,j);
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // reset state at split
            in_f = 1;
            Reset(X_split, t_split);
            Y = Y_split;
            Z = Z_split;
            xi = 0.0;
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // continue until exit
            do {
              counter += D;
              Increment_Update();
              VR_CV_Step();
              //N_steps+=1/split;
              #ifdef PLOT
              if(np < 10){
                file_plot = fopen(filename, "a");
                for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
                fprintf(file_plot,"\n");
                fclose(file_plot);
              }
              #endif
              if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
            } while (in_f);
            if(status == time_out){
              summ_nvr += (Y*bvp.p.Value(X,0.0)+Z);
              sqsumm_nvr += pow(Y*bvp.p.Value(X,0.0)+Z,2);
              summ_vr += (Y*bvp.p.Value(X,0.0)+Z+xi);
              sqsumm_vr += pow(Y*bvp.p.Value(X,0.0)+Z+xi,2);
            } else {
              summ_nvr += (Y*bvp.g.Value(E_P,t)+Z);
              sqsumm_nvr += pow(Y*bvp.p.Value(E_P,t)+Z,2);
              summ_vr += (Y*bvp.g.Value(E_P,t)+Z+xi);
              sqsumm_vr += pow(Y*bvp.p.Value(E_P,t)+Z+xi,2);
            }
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
          }
          var_vr = pow(summ_vr/10.0,2.0) - sqsumm_vr/10.0;
          var_nvr = pow(summ_nvr/10.0,2.0) - sqsumm_nvr/10.0;
          split = std::min((int)((1<<l)*var_vr/var_nvr), 1<<l);
          split_f = (double) split;
          for (int s=0; s<split-10; s++) {
            #ifdef PLOT
            if(np < 10){
              sprintf(filename,"Output/trayectories/split_%d.csv",s);
              file_plot = fopen(filename, "w");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"X_%u,E_P_%u,N_%u,",j,j,j);
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // reset state at split
            in_f = 1;
            Reset(X_split, t_split);
            Y = Y_split;
            Z = Z_split;
            xi = 0.0;
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
	          // continue until exit
            do {
              counter += D;
              Increment_Update();
              VR_CV_Step();
              //N_steps+=1/split;
              #ifdef PLOT
              if(np < 10){
                file_plot = fopen(filename, "a");
                for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
                fprintf(file_plot,"\n");
                fclose(file_plot);
              }
              #endif
              if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
            } while (in_f);
            if(status == time_out){
              summ_vr += (Y*bvp.p.Value(X,0.0)+Z+xi);
              sqsumm_vr += pow(Y*bvp.p.Value(X,0.0)+Z+xi,2);
            } else {
              summ_vr += (Y*bvp.g.Value(E_P,t)+Z+xi);
              sqsumm_vr += pow(Y*bvp.p.Value(E_P,t)+Z+xi,2);
            }
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",X(j),E_P(j),N(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
          }
          Pf = summ_vr/split;
          var_vr = pow(summ_vr/split,2.0) - sqsumm_vr/split;
        }
      }

      if (in_c) {
        t_split = tc;
        X_split = Xc;
        Y_split = Yc;
        Z_split = Zc;
        xi_split = 0.0;
        split_c = (double) split;
        for (int s=0; s<split; s++) {
	        // reset state at split
          #ifdef PLOT
          if(np < 10){
            sprintf(filename,"Output/trayectories/split_%d.csv",s);
            file_plot = fopen(filename, "w");
            for(j = 0; j < X0.size(); j++) fprintf(file_plot,"X_%u,E_P_%u,N_%u,",j,j,j);
            fprintf(file_plot,"\n");
            fclose(file_plot);
          }
          #endif
          in_c = 1;
          Resetc(X_split, t_split);
          Yc = Y_split;
          Zc = Z_split;
          xic = xi_split;
          #ifdef PLOT
          if(np < 10){
            file_plot = fopen(filename, "a");
            for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
            fprintf(file_plot,"\n");
            fclose(file_plot);
          }
          #endif
	        // continue until exit
          do {
            counter += D;
            Incrementc_Update();

            Stepc();
            #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
            #endif
            if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
          } while (in_c);
          if(status == time_out){
                Pc += (Yc*bvp.p.Value(Xc,0.0)+Zc+xic)/ ((double) split);
                sqsumm_vr += pow(Yc*bvp.p.Value(Xc,0.0)+Zc+xic,2)/ 
                ((double) split);
          } else {
                Xc= E_Pc;
                Pc += (Yc*bvp.g.Value(E_Pc,tc)+Zc+xic)/ ((double) split);
                sqsumm_vr += pow(Yc*bvp.p.Value(Xc,tc)+Zc+xic,2)/ 
                ((double) split);
          }
          #ifdef PLOT
            if(np < 10){
              file_plot = fopen(filename, "a");
              for(j = 0; j < X0.size(); j++) fprintf(file_plot,"%f,%f,%f,",Xc(j),E_Pc(j),Nc(j));
              fprintf(file_plot,"\n");
              fclose(file_plot);
            }
          #endif
        }
      }

    }
    //printf("Split c %f Split f %f \n", split_c, split_f);
    dP = Pf-Pc;
    sol = sol_0 + dP;
    sums[0] += counter; // add number of RNG calls
    sums[1] += dP;
    sums[2] += dP*dP;
    sums[3] += dP*dP*dP;
    sums[4] += dP*dP*dP*dP;
    sums[5] += split_f;
    sums[6] += split_c;
    sums[7] += pow(sol-sol_a,2.0); //MSE
    sums[8] += N_steps;
    if(split_f > 0.1) sums[9] += var_vr;
    #ifdef PLOT
        if(np  < 10){
          sprintf(filename,"python3 trayectoriesmlvl2D.py %d",(int)(split_f+split_c));
          system(filename);
        }
    #endif
  }
}
*/
void MLVLSolver::Solve_l(Eigen::VectorXd X0, double t_start, double noise_mean, double noise_var){
  int16_t split, in_f, in_c, set_f, set_c;
  uint32_t counter = 0;
  double t_split, Y_split, Z_split, xi_split, split_c, split_f, N_steps;
  Eigen::VectorXd X_split;
  increment.resize(X0.size());
  incrementc.resize(X0.size());
  N.resize(X0.size());
  E_P.resize(X0.size());
  Nc.resize(X0.size());
  E_Pc.resize(X0.size());
  D = X0.size();
  h = h0 / pow(M,l);
  sqrth = sqrt(h);
  hc = h0 / pow(M,(l-1));
  sqrthc = sqrt(hc);
  //printf(" **** l, Nl = %u, %u **** h = %f hc = %f\n",l,dNl[l],h,hc);
  sol_a = bvp.u.Value(X0,t_start);
  sol_0 = 0.0;
  if (l == 0){
      sol_0 = 0;
  } else {
      for(int itl = 0; itl < l-1; itl++) {
        sol_0 += ml[itl];
      }
  }
  option = 1;
  for (unsigned int k=0; k<10; k++) sums[k] = 0.0;
  if (option==1 || option==3) {
    //split = 4 * (1<<l);
    split = 1<<l;
    //split = 0;
  }
  else if (option==2) {
    split = 1;
  }
  
  for (unsigned int np = 0; np < dNl[l]; np++) {
    in_f = 1;
    in_c = 1;
    set_f = 0;
    set_c = 0;
    Pf = 0.0;
    Pc = 0.0;
    counter = 0;
    split_c = 0.0;
    split_f = 0.0;
    N_steps = 0;
    // level 0
    Reset(X0, t_start);
    Resetc(X0, t_start);
    if (l==0) {
      do {
        Increment_Update();
        counter += D;
        VR_CV_Step();
        N_steps ++;
        if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
      } while (in_f);
      if(status == time_out){
          Pf = Y*bvp.p.Value(X,0.0)+Z+xi;
      } else {
          #ifdef NOISYBOUNDARY
          Pf = Y*(bvp.g.Value(E_P,t) + noise_mean + 
          gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Z+xi;
          #else
          Pf = Y*bvp.g.Value(E_P,t)+Z+xi;
          #endif
      }
    }

      // level l>0

    else {

      Reset(X0, t_start);
      Resetc(X0, t_start);
      do {

        for (int d=0; d<D; d++) incrementc(d) = 0.0;
        for (int m=0; m<M; m++){
          counter += D;
          Increment_Update();
          incrementc += increment;
          VR_CV_Step();
          N_steps ++;
          if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
          if ( (!in_f) && (!set_f) ) {
            set_f = 1;
            if(status == time_out){
                Pf = Y*bvp.p.Value(X,0.0)+Z+xi;
            } else {
                #ifdef NOISYBOUNDARY
                Pf = Y*(bvp.g.Value(E_P,t) + noise_mean + 
                gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Z+xi;
                #else
                Pf = Y*bvp.g.Value(E_P,t)+Z+xi;
                #endif
            }

            break;
          }
      	}
        Stepc();
        if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
        if ( (!in_c) && (!set_c) ) {
          set_c = 1;
          if(status == time_out){
                Pc = Yc*bvp.p.Value(Xc,0.0)+Zc+xic;
          } else {
                #ifdef NOISYBOUNDARY
                Pc = Yc*(bvp.g.Value(E_Pc,tc)+noise_mean + 
                gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Zc+xic;
                #else
                Pc = Yc*bvp.g.Value(E_Pc,tc)+Zc+xic;
                #endif
          }

        }
      } while (in_f && in_c);
      // split continuation paths
      #ifdef NOISYPSEUDOSPLIT
      if (in_f) {
        Pf = Y*(bvp.u.Value(E_P,t) + noise_mean + 
          gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Z+xi;
      }
      if (in_c) {
        Pc = Yc*(bvp.u.Value(E_Pc,tc) + noise_mean + 
          gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Zc+xic;
      }
      #else
      if (in_f) {
        t_split = t;
        X_split = X;
        Y_split = Y;
        Z_split = Z;
        xi_split = xi;
        split_f = (double) split;
        for (int s=0; s<split; s++) {
	        // reset state at split
          in_f = 1;
          Reset(X_split, t_split);
          Y = Y_split;
          Z = Z_split;
          xi = xi_split;
	        // continue until exit
          do {
            counter += D;
            Increment_Update();
            VR_CV_Step();
            N_steps+=1/split;
            if(!Inside(X,N,E_P,t,ji_t,sqrth)) in_f = 0;
          } while (in_f);
          if(status == time_out){
                Pf += (Y*bvp.p.Value(X,0.0)+Z+xi)/((double) split);
          } else {
                #ifdef NOISYBOUNDARY
                Pf += (Y*(bvp.g.Value(E_P,t) + noise_mean + 
                gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Z+xi)/((double) split);
                #else
                Pf += (Y*bvp.g.Value(X,t)+Z+xi)/((double) split);
                #endif
          }

        }
      }

      if (in_c) {
        t_split = tc;
        X_split = Xc;
        Y_split = Yc;
        Z_split = Zc;
        xi_split = xic;
        split_c = (double) split;
        for (int s=0; s<split; s++) {
	        // reset state at split
          in_c = 1;
          Resetc(X_split, t_split);
          Yc = Y_split;
          Zc = Z_split;
          xic = xi_split;
	        // continue until exit
          do {
            counter += D;
            Incrementc_Update();

            Stepc();
            if(!Inside(Xc,Nc,E_Pc,tc,ji_tc,sqrthc)) in_c = 0;
          } while (in_c);
          if(status == time_out){
                Pc += (Yc*bvp.p.Value(Xc,0.0)+Zc+xic)/ ((double) split);
          } else {
                #ifdef NOISYBOUNDARY
                Pc += (Yc*(bvp.g.Value(E_Pc,tc)+noise_mean + 
                gsl_ran_gaussian_ziggurat(rng,1)*sqrt(noise_var))+Zc+xic)/
                ((double) split);
                #else
                Pc += (Yc*(bvp.g.Value(E_Pc,tc)+Zc+xic))/
                ((double) split);
                #endif
          }
        }
      }
      #endif

    }
    //printf("Split c %f Split f %f \n", split_c, split_f);
    dP = Pf-Pc;
    sol = sol_0 + dP;
    sums[0] += counter; // add number of RNG calls
    sums[1] += dP;
    sums[2] += dP*dP;
    sums[3] += dP*dP*dP;
    sums[4] += dP*dP*dP*dP;
    sums[5] += split_f;
    sums[6] += split_c;
    sums[7] += pow(sol-sol_a,2.0); //MSE
    sums[8] += N_steps; 
  }
}

void MLVLSolver::Regression(unsigned int N_sample, double *x, double *y, double &a, double &b){
  double sum0=0.0, sum1=0.0, sum2=0.0, sumy0=0.0, sumy1=0.0;

  for (uint32_t i=0; i<N_sample; i++) {
    sum0  += 1.0;
    sum1  += x[i];
    sum2  += x[i]*x[i];

    sumy0 += y[i];
    sumy1 += y[i]*x[i];
  }

  a = (sum0*sumy1 - sum1*sumy0) / (sum0*sum2 - sum1*sum1);
  b = (sum2*sumy0 - sum1*sumy1) / (sum0*sum2 - sum1*sum1);
}
double MLVLSolver::Solve(Eigen::VectorXd X0, int Lmin, int Lmax, int N0, double eps,
                           double alpha_0,double beta_0,double gamma_0){
      
      return Solve(X0, INFINITY, Lmin, Lmax, N0, eps, alpha_0,beta_0,gamma_0);
}

double MLVLSolver::Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
  double alpha_0,double beta_0,double gamma_0){
  double suml[6][21], rem, P;
  double x[21], y[21], alpha, beta, gamma, sum, theta;
  int    converged;
  int    diag = 0;  // diagnostics, set to 0 for none 
  //
  // check input parameters
  //
  E_P.resize(X0.size());
  N.resize(X0.size());
  if (bvp.boundary.Dist(params,X0,E_P,N) > -1e-06) return bvp.g.Value(X0,t_start);
  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0,alpha_0);
  beta  = fmax(0.0,beta_0);
  gamma = fmax(0.0,gamma_0);
  theta = 0.25;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0,(double)l*gamma);
    NlCl[l] = 0.0;
    dNl[l]  = 0;
    for(int n=0; n<6; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //
    if(diag) printf("Computing L = %d \n With N ", L);
    for (l=0; l<=L; l++) {
      if (diag) printf("N_%d %u ",l,dNl[l]);
      if (dNl[l]>0) {
        Solve_l(X0, t_start);
        suml[0][l] += (double) dNl[l]; //Number of trayectories
        suml[1][l] += sums[1]; //Correction
        for (int itl=0; itl<=l; itl++) ml[itl] = suml[1][itl]/suml[0][itl];
        suml[2][l] += sums[2]; //Square of correction
        suml[3][l] += sums[7]; //MSE
        suml[4][l] += sums[5]; //Fine Splitted trayectories
        suml[5][l] += sums[6]; //Coarsed Splitted trayectories
        NlCl[l]    += sums[0]; //sum total cost
      }
    }
    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0;

    for (int l=0; l<=L; l++) {
      ml[l] = suml[1][l]/suml[0][l];//Average
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
      if (gamma_0 <= 0.0) Cl[l] = NlCl[l] / suml[0][l];//Average cost

      if (l>1) {
        ml[l] = suml[1][l]/suml[0][l];//Average
        //Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
        //ml[l] = fmaxf(ml[l],  0.5*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5*Vl[l-1]/powf(2.0,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }
    
    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0-theta)*eps*eps)
                     - suml[0][l]));
      
    }
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(fabs(ml[l]));
        //if (diag) printf("m[%d] = %f\n",l,ml[l]);
      }
      Regression(L,x,y,alpha,sum);
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      Regression(L,x,y,beta,sum);
      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      Regression(L,x,y,gamma,sum);
      if (diag) printf(" gamma = %f \n",gamma);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
    for (int l=0; l<=L; l++) sum += fmaxf(0.0, dNl[l]-0.01*suml[0][l]);
    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      rem = fabs(ml[L]) / (powf(2.0f,alpha)-1.0f);
      //if (diag) std::cout << rem << "alpha is " << alpha << std::endl;
      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf(fmaxf( 0.0f, 
                          sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
  
    P = 0.0;
    var = 0.0;
    for (int l=0; l<=L; l++) {
      P    += suml[1][l]/suml[0][l];
      var += Vl[l]/suml[0][l];
      Nl[l] = suml[0][l];
      Cl[l] = NlCl[l] / Nl[l];
      mse[l] = suml[3][l]/suml[0][l];
    }
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  return P;
}
double MLVLSolver::Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
  double alpha_0,double beta_0,double gamma_0,std::string fn_results, std::string fn_nl){
  double suml[9][21], rem, P;
  double x[21], y[21], alpha, beta, gamma, sum, theta;
  int    converged;
  int    diag = 0;  // diagnostics, set to 0 for none
  unsigned int t0,t1;
  //
  // check input parameters
  //
  FILE *fp;
  fp = fopen(fn_results.c_str(),"w");
  fprintf(fp,"L,sol_a,sol_n,err,rerr,mse,var\n");
  fclose(fp);
  fp = fopen(fn_nl.c_str(),"w");
  fprintf(fp,"l,hf,hc,ml,Vl,Cl,Nl,Nf,Nc,APL,var_split,Ctime\n");
  fclose(fp);
  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0,alpha_0);
  beta  = fmax(0.0,beta_0);
  gamma = fmax(0.0,gamma_0);
  theta = 0.25;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0,(double)l*gamma);
    NlCl[l] = 0.0;
    dNl[l]  = 0;
    for(int n=0; n<9; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //
    if(diag) printf("Computing L = %d \n With N ", L);
    for (l=0; l<=L; l++) {
      //if (diag) printf("N_%d %u ",l,dNl[l]);
      if (dNl[l]>0) {
        t0 = clock();
        Solve_l(X0, t_start);
        t1 = clock();
        suml[0][l] += (double) dNl[l]; //Number of trayectories
        suml[1][l] += sums[1]; //Correction
        for (int itl=0; itl<=l; itl++) ml[itl] = suml[1][itl]/suml[0][itl];
        suml[2][l] += sums[2]; //Square of correction
        for (int itl=0; itl<=l; itl++) Vl[itl] = suml[2][l]/suml[0][l] - ml[l]*ml[l];
        suml[3][l] += sums[7]; //MSE
        suml[4][l] += sums[5]; //Fine Splitted trayectories
        suml[5][l] += sums[6]; //Coarsed Splitted trayectories
        suml[6][l] += sums[8]; //Total number of steps
        NlCl[l]    += sums[0]; //sum total cost
        fp = fopen(fn_nl.c_str(),"a");
        fprintf(fp,"%u,%f,%f,%f,%f,%f,%u,%f,%f,%f,%f,%u\n",
                l, h0/pow(M,l),h0/pow(M,l-1), ml[l], Vl[l], 
                NlCl[l],dNl[l],suml[4][l],suml[5][l],suml[6][l]/dNl[l]
                ,suml[7][l],t1-t0);
        fclose(fp);
      }
    }
    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0;

    for (int l=0; l<=L; l++) {
      ml[l] = suml[1][l]/suml[0][l];//Average
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
      if (gamma_0 <= 0.0) Cl[l] = NlCl[l] / suml[0][l];//Average cost

      if (l>1) {
        ml[l] = suml[1][l]/suml[0][l];//Average
        //Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
        //ml[l] = fmaxf(ml[l],  0.5*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5*Vl[l-1]/powf(2.0f,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }
    
    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0-theta)*eps*eps)
                     - suml[0][l]));
      
    }
    fprintf(fp,"\n");
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(fabs(ml[l]));
      }
      Regression(L,x,y,alpha,sum);
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      Regression(L,x,y,beta,sum);
      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      Regression(L,x,y,gamma,sum);
      if (diag) printf(" gamma = %f \n",gamma);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
    for (int l=0; l<=L; l++) sum += fmaxf(0.0, dNl[l]-0.01*suml[0][l]);
    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      rem = fabs(ml[L]) / (powf(2.0f,alpha)-1.0f);
      //std::cout << rem << "alpha is " << alpha << std::endl;
      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf(fmaxf( 0.0f, 
                          sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
    P = 0.0;
    var = 0.0;
    for (int l=0; l<=L; l++) {
      P    += suml[1][l]/suml[0][l];
      var += Vl[l];
      Nl[l] = suml[0][l];
      Cl[l] = NlCl[l] / Nl[l];
      mse[l] = suml[3][l]/suml[0][l];
    }
    fp = fopen(fn_results.c_str(),"a");
    fprintf(fp,"%u,%e,%e,%e,%e,%e,%e\n", L, sol_a,
            P,fabs(sol_a-P), fabs(sol_a-P)/sol_a, mse[L],var);
    fclose(fp);
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  return P;
}
double MLVLSolver::Solve(Eigen::VectorXd X0, double t_start, int Lmin, int Lmax, int N0, double eps,
  double noise_mean, double noise_variance, double alpha_0,double beta_0,double gamma_0,
  std::string fn_results, std::string fn_nl){
  double suml[8][21], rem, P;
  double x[21], y[21], alpha, beta, gamma, sum, theta;
  int    converged;
  int    diag = 0;  // diagnostics, set to 0 for none
  unsigned int t0,t1;
  //
  // check input parameters
  //
  FILE *fp;
  fp = fopen(fn_results.c_str(),"w");
  fprintf(fp,"L,sol_a,sol_n,err,rerr,mse,var,n_var,n_mean\n");
  fclose(fp);
  fp = fopen(fn_nl.c_str(),"w");
  fprintf(fp,"l,hf,hc,ml,Vl,Cl,Nl,Nf,Nc,APL,Ctime,n_var,n_mean\n");
  fclose(fp);
  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0,alpha_0);
  beta  = fmax(0.0,beta_0);
  gamma = fmax(0.0,gamma_0);
  theta = 0.25;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0,(double)l*gamma);
    NlCl[l] = 0.0;
    dNl[l]  = 0;
    for(int n=0; n<8; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //
    if(diag) printf("Computing L = %d \n With N ", L);
    for (l=0; l<=L; l++) {
      //if (diag) printf("N_%d %u ",l,dNl[l]);
      if (dNl[l]>0) {
        t0 = clock();
        Solve_l(X0, t_start, noise_mean, noise_variance);
        t1 = clock();
        suml[0][l] += (double) dNl[l]; //Number of trayectories
        suml[1][l] += sums[1]; //Correction
        for (int itl=0; itl<=l; itl++) ml[itl] = suml[1][itl]/suml[0][itl];
        suml[2][l] += sums[2]; //Square of correction
        for (int itl=0; itl<=l; itl++) Vl[itl] = suml[2][l]/suml[0][l] - ml[l]*ml[l];
        suml[3][l] += sums[7]; //MSE
        suml[4][l] += sums[5]; //Fine Splitted trayectories
        suml[5][l] += sums[6]; //Coarsed Splitted trayectories
        suml[6][l] += sums[8]; //Total number of steps
        NlCl[l]    += sums[0]; //sum total cost
        fp = fopen(fn_nl.c_str(),"a");
        fprintf(fp,"%u,%f,%f,%f,%f,%f,%u,%f,%f,%f,%u,%f,%f\n",
                l, h0/pow(M,l),h0/pow(M,l-1), ml[l], Vl[l], 
                NlCl[l],dNl[l],suml[4][l],suml[5][l],suml[6][l]/dNl[l]
                ,t1-t0,noise_variance,noise_mean);
        fclose(fp);
      }
    }
    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0;

    for (int l=0; l<=L; l++) {
      ml[l] = suml[1][l]/suml[0][l];//Average
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
      if (gamma_0 <= 0.0) Cl[l] = NlCl[l] / suml[0][l];//Average cost

      if (l>1) {
        ml[l] = suml[1][l]/suml[0][l];//Average
        //Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0);//variance
        //ml[l] = fmaxf(ml[l],  0.5*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5*Vl[l-1]/powf(2.0f,beta));
      }

      sum += sqrtf(Vl[l]*Cl[l]);
    }
    
    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0-theta)*eps*eps)
                     - suml[0][l]));
      
    }
    fprintf(fp,"\n");
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(fabs(ml[l]));
      }
      Regression(L,x,y,alpha,sum);
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      Regression(L,x,y,beta,sum);
      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      Regression(L,x,y,gamma,sum);
      if (diag) printf(" gamma = %f \n",gamma);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
    for (int l=0; l<=L; l++) sum += fmaxf(0.0, dNl[l]-0.01*suml[0][l]);
    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      rem = fabs(ml[L]) / (powf(2.0f,alpha)-1.0f);
      //std::cout << rem << "alpha is " << alpha << std::endl;
      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf(fmaxf( 0.0f, 
                          sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
    P = 0.0;
    var = 0.0;
    sol_a = bvp.u.Value(X0,t_start);
    for (int l=0; l<=L; l++) {
      P    += suml[1][l]/suml[0][l];
      var += Vl[l];
      Nl[l] = suml[0][l];
      Cl[l] = NlCl[l] / Nl[l];
      mse[l] = suml[3][l]/suml[0][l];
    }
    fp = fopen(fn_results.c_str(),"a");
    fprintf(fp,"%u,%e,%e,%e,%e,%e,%e,%e,%e\n", L, sol_a,
            P,fabs(sol_a-P), fabs(sol_a-P)/sol_a, mse[L],var,
            noise_variance, noise_mean);
    fclose(fp);
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  return P;
}
void MLVLSolver::Test(Eigen::VectorXd X0, uint32_t N_test, uint32_t L_test, uint32_t N0, 
      double *Eps, uint32_t Lmin, uint32_t Lmax, std::string fname){
    Test(X0,INFINITY,N_test, L_test, N0, Eps, Lmin, Lmax, fname);
}
 void MLVLSolver::Test(Eigen::VectorXd X0, double t_start, uint32_t N_test, uint32_t L_test,
     uint32_t N0, double *Eps, uint32_t Lmin, uint32_t Lmax, std::string fname){
  //
  // first, convergence tests
  //

  // current date/time based on current system
  time_t now = time(NULL);
  char *date = ctime(&now);
  int len = strlen(date);
  date[len-1] = ' ';

  printf("\n");
  printf("**********************************************************\n");
  printf("*** MLMC file version 0.9     produced by              ***\n");
  printf("*** C++ mlmc_test on %s         ***\n",date);
  printf("**********************************************************\n");
  printf("\n");
  printf("**********************************************************\n");
  printf("*** Convergence tests, kurtosis, telescoping sum check ***\n");
  printf("*** using N =%7d samples                           ***\n",N_test);
  printf("**********************************************************\n");
  printf("\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  printf("    kurtosis     check        cost \n--------------------------");
  printf("-------------------------------------------------------------\n");
  
  float *cost = (float *)malloc((L_test+1)*sizeof(float));
  float *del1 = (float *)malloc((L_test+1)*sizeof(float));
  float *del2 = (float *)malloc((L_test+1)*sizeof(float));
  float *var1 = (float *)malloc((L_test+1)*sizeof(float));
  float *var2 = (float *)malloc((L_test+1)*sizeof(float));
  float *chk1 = (float *)malloc((L_test+1)*sizeof(float));
  float *kur1 = (float *)malloc((L_test+1)*sizeof(float));

  for (l=0; l<=L_test; l++) {
    dNl[l] = N_test;
    Solve_l(X0, t_start);

    for (int m=0; m<7; m++) sums[m] = sums[m]/N_test;

    if (M>0) cost[l] = powf((float)M,(float)l);
    else cost[l] = sums[0];
    del1[l] = sums[1];
    del2[l] = sums[5];
    var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
    var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

    kur1[l]  = (      sums[4]
                - 4.0*sums[3]*sums[1]
                + 6.0*sums[2]*sums[1]*sums[1]
                - 3.0*sums[1]*sums[1]*sums[1]*sums[1] )
             / (var1[l]*var1[l]);

    if (l==0)
      chk1[l] = 0.0f;
    else
      chk1[l] = sqrtf((float) N_test) * 
                fabsf(  del1[l]  +       del2[l-1]  -       del2[l] )
         / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

    printf("%2d  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",
    l,del1[l],del2[l],var1[l],var2[l],kur1[l],chk1[l],cost[l]/N_test);
    L = l;
  }

//
// print out a warning if kurtosis or consistency check looks bad
//

  if (kur1[L] > 100.0f) {
    printf("\n WARNING: kurtosis on finest level = %f \n",kur1[L]);
    printf(" indicates MLMC correction dominated by a few rare paths; \n");
    printf(" for information on the connection to variance of sample variances,\n");
    printf(" see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  float max_chk = 0.0f;
  for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
  if (max_chk > 1.0f) {
    printf("\n WARNING: maximum consistency error = %f \n",max_chk);
    printf(" indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

//
// use linear regression to estimate alpha, beta, gamma
//

  double alpha, beta, gamma, foo;
  double *x, *y;
  x = new double[L];
  y = new double[L];

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(fabsf(del1[l]));
  } 
  Regression(L,x,y,alpha,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(var1[l]);
  } 
  Regression(L,x,y,beta,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = log2f(cost[l]);
  } 
  Regression(L,x,y,gamma,foo);

  printf("\n******************************************************\n");
  printf("*** Linear regression estimates of MLMC parameters ***\n");
  printf("******************************************************\n");
  printf("\n alpha = %f  (exponent for MLMC weak convergence)\n",alpha);
  printf(" beta  = %f  (exponent for MLMC variance) \n",beta);
  printf(" gamma = %f  (exponent for MLMC cost) \n",gamma);

//
// second, mlmc complexity tests
//

  printf("\n");
  printf("***************************** \n");
  printf("*** MLMC complexity tests *** \n");
  printf("***************************** \n\n");
  printf("  eps       value   a.error  r.error   mlmc_cost   std_cost  savings     N_l \n");
  printf("--------------------------------------------------------- \n");
 
  int i=0;

  while (Eps[i]>0) {
    double eps = Eps[i++];

    double P = Solve(X0, t_start, Lmin, Lmax, N0, eps, alpha, beta, gamma);

    double std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

    for (uint32_t l=0; l<=Lmax; l++) {
      if (Nl[l]>0) {
        // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
        mlmc_cost += Nl[l]*Cl[l];
        if (l<=L)
          std_cost = var2[l]*Cl[l] / ((1.0f-theta)*eps*eps);
        else
          std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
      }
    }

    printf("%.4f  %.4e %.3e %.3e %.3e  %.3e  %7.2f ",
	          eps, P, P-bvp.u.Value(X0,t_start), fabs(P-bvp.u.Value(X0,t_start))/bvp.u.Value(X0, t_start), mlmc_cost, std_cost, std_cost/mlmc_cost);
    for (int l=0; Nl[l]>0; l++) printf("%9d",Nl[l]);
    printf("\n");
  }
  printf("\n");
  delete x;
  delete y;
}