#include "stencil.hpp"

Stencil::Stencil(void){
    //All variables of the class are emptied
    pos_north.clear(); index_north.clear();
    pos_south.clear(); index_south.clear();
    pos_east.clear(); index_east.clear();
    pos_west.clear(); index_west.clear();
}

void Stencil::Init(std::map<direction, std::vector<int> > s_index,
                   std::map<direction, std::vector<double> > s_x,
                   std::map<direction, std::vector<double> > s_y,
                   double *parameters){
    pos_north.clear(); index_north.clear(); G_north.clear();
    pos_south.clear(); index_south.clear(); G_south.clear();
    pos_east.clear(); index_east.clear(); G_east.clear();
    pos_west.clear(); index_west.clear();G_west.clear();
    ipsi_north.resize(0,0); ipsi_south.resize(0,0);
    ipsi_east.resize(0,0); ipsi_west.resize(0,0);
    std::map<int,Eigen::VectorXd> sten_map;
    Eigen::VectorXd vaux;
    vaux.resize(2);
    for(int i = 0; i < 4; i++) stencil_parameters[i] = parameters[i];
    if(s_index[North].size() > 0){
        for(int i = 0; i < (int)s_index[North].size(); i++){
            vaux[0] = s_x[North][i];
            vaux[1] = s_y[North][i];
            sten_map[s_index[North][i]] = vaux;
            index_north.push_back(s_index[North][i]);
        }
        std::vector<int>::iterator v_north = std::unique(index_north.begin(), index_north.end());
        index_north.erase(v_north, index_north.end());
        for(unsigned int i = 0; i < index_north.size(); i++) pos_north.push_back(sten_map[index_north[i]]);
        for(unsigned int i = 0; i < pos_north.size(); i++) G_north.push_back(0.0);
    }
    sten_map.clear();

    if(s_index[South].size() > 0){
        for(int i = 0; i < (int)s_index[South].size(); i++){
            vaux[0] = s_x[South][i];
            vaux[1] = s_y[South][i];
            sten_map[s_index[South][i]] = vaux;
            index_south.push_back(s_index[South][i]);
        }
        std::vector<int>::iterator v_south = std::unique(index_south.begin(), index_south.end());
        index_south.erase(v_south, index_south.end());
        for(unsigned int i = 0; i < index_south.size(); i++) pos_south.push_back(sten_map[index_south[i]]);
        for(unsigned int i = 0; i < pos_south.size(); i++) G_south.push_back(0.0);
    }
    sten_map.clear();

    if(s_index[East].size() > 0){
        for(int i = 0; i < (int)s_index[East].size(); i++){
            vaux[0] = s_x[East][i];
            vaux[1] = s_y[East][i];
            sten_map[s_index[East][i]] = vaux;
            index_east.push_back(s_index[East][i]);
        }
        std::vector<int>::iterator v_east = std::unique(index_east.begin(), index_east.end());
        index_east.erase(v_east, index_east.end());
        for(unsigned int i = 0; i < index_east.size(); i++) pos_east.push_back(sten_map[index_east[i]]);
        for(unsigned int i = 0; i < pos_east.size(); i++) G_east.push_back(0.0);
    }
    sten_map.clear();

    if(s_index[West].size() > 0){
        for(int i = 0; i < (int)s_index[West].size(); i++){
            vaux[0] = s_x[West][i];
            vaux[1] = s_y[West][i];
            sten_map[s_index[West][i]] = vaux;
            index_west.push_back(s_index[West][i]);
        }
        std::vector<int>::iterator v_west = std::unique(index_west.begin(), index_west.end());
        index_west.erase(v_west, index_west.end()); 
        for(unsigned int i = 0; i < index_west.size(); i++) pos_west.push_back(sten_map[index_west[i]]);
        for(unsigned int i = 0; i < pos_west.size(); i++) G_west.push_back(0.0);
    }
    sten_map.clear();

    counter_north = 0; counter_south = 0;
    counter_east = 0; counter_west = 0;
}
void Stencil::Reset(void){
    G_north.clear();
    for(unsigned int i = 0; i < pos_north.size(); i++) G_north.push_back(0.0);
    G_south.clear();
    for(unsigned int i = 0; i < pos_south.size(); i++) G_south.push_back(0.0);
    G_east.clear();
    for(unsigned int i = 0; i < pos_east.size(); i++) G_east.push_back(0.0);
    G_west.clear();
    for(unsigned int i = 0; i < pos_west.size(); i++) G_west.push_back(0.0);
    counter_north = 0; counter_south = 0;
    counter_east = 0; counter_west = 0;
}

Eigen::MatrixXd Stencil::Compute_ipsi(std::vector<Eigen::VectorXd> & sten_position,
                  BVP bvp,
                  double c2){
    //Auxiliary vectors 
    Eigen::VectorXd vaux, sincrement;
    //Psi Matrix, its inverse and Identity are created 
    Eigen::MatrixXd Psi, iPsi, I;
    //Matrix are resized
    Psi.resize(sten_position.size(), sten_position.size());
    I.resize(sten_position.size(), sten_position.size());
    I.setIdentity();
    //And fullfilled
    for(unsigned int i = 0; i < sten_position.size(); i ++){
        for(unsigned int j = 0; j < sten_position.size(); j ++){
            Psi(i,j) = bvp.rbf.Value(sten_position[i], sten_position[j], c2);
        }
    }
    double err, cond;
    iPsi= Psi.inverse();
    Eigen::BiCGSTAB<Eigen::MatrixXd> CGS;
    CGS.compute(Psi);
    iPsi = CGS.solveWithGuess(I,iPsi);
    err = CGS.error();
    cond = Psi.norm()*iPsi.norm();
    //std::cout << "Inverse matrix computed with error "<< err <<" and condition number" << cond << std::endl;
    Psi.resize(0,0);
    I.resize(0,0);
    return iPsi;
}

bool Stencil::AreSame(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

void Stencil::Compute_ipsi(BVP bvp, double c2){
	if(pos_north.size() > 0){
		ipsi_north = Compute_ipsi(pos_north, bvp, c2);
	}
	if(pos_south.size() > 0){
		ipsi_south = Compute_ipsi(pos_south, bvp, c2);
	}
	if(pos_east.size() > 0){
		ipsi_east = Compute_ipsi(pos_east, bvp, c2);
	}
	if(pos_west.size() > 0){
		ipsi_west = Compute_ipsi(pos_west, bvp, c2);
	}
}

int Stencil::G_update(Eigen::VectorXd X, double Y, BVP bvp, double c2){
	double H_ij;
	if(AreSame(X(1),stencil_parameters[1])){
		//South stencil
        //std::cout << "\t south \n";
		counter_south ++;
		for(unsigned int i = 0; i < G_south.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_south.size(); j++){
				H_ij += ipsi_south(j,i)*bvp.rbf.Value(pos_south[j],X,c2);
			}
			G_south[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(1),stencil_parameters[3])){
		//North stencil
        //std::cout << "\t north \n";
		counter_north ++;
		for(unsigned int i = 0; i < G_north.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_north.size(); j++){
				H_ij += ipsi_north(j,i)*bvp.rbf.Value(pos_north[j],X,c2);
			}
			G_north[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(0),stencil_parameters[0])){
		//West stencil
        //std::cout << "\t west \n";
		counter_west ++;
		for(unsigned int i = 0; i < G_west.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_west.size(); j++){
				H_ij += ipsi_west(j,i)*bvp.rbf.Value(pos_west[j],X,c2);
			}
			G_west[i] += -H_ij*Y;
		}
		return 1;
	}
	if(AreSame(X(0),stencil_parameters[2])){
		//East stencil
        //std::cout << "\t east \n";
		counter_east ++;
		for(unsigned int i = 0; i < G_east.size(); i++){
            H_ij = 0.0f;
			for(unsigned int j = 0; j < G_east.size(); j++){
				H_ij += ipsi_east(j,i)*bvp.rbf.Value(pos_east[j],X,c2);
			}
			G_east[i] += -H_ij*Y;
		}
		return 1;
	}
	std::cout << "Something went wrong computing H.";
	return 0;
}

double Stencil::Test_Interpolator(Eigen::VectorXd X, BVP bvp, double c2){
    double uH_ij;
    if(AreSame(X(1),stencil_parameters[1])){
        //South stencil
        //std::cout << "\t south \n";
        counter_south ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_south.size(); i++){
            for(unsigned int j = 0; j < G_south.size(); j++){
                uH_ij += ipsi_south(j,i)*bvp.rbf.Value(pos_south[j],X,c2)*bvp.u.Value(pos_south[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(1),stencil_parameters[3])){
        //North stencil
        //std::cout << "\t north \n";
        counter_north ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_north.size(); i++){
            for(unsigned int j = 0; j < G_north.size(); j++){
                uH_ij += ipsi_north(j,i)*bvp.rbf.Value(pos_north[j],X,c2)*bvp.u.Value(pos_north[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(0),stencil_parameters[0])){
        //West stencil
        //std::cout << "\t west \n";
        counter_west ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_west.size(); i++){
            for(unsigned int j = 0; j < G_west.size(); j++){
                uH_ij += ipsi_west(j,i)*bvp.rbf.Value(pos_west[j],X,c2)*bvp.u.Value(pos_west[i],0.0f);
            }
        }
        return uH_ij;
    }
    if(AreSame(X(0),stencil_parameters[2])){
        //East stencil
        //std::cout << "\t east \n";
        counter_east ++;
        uH_ij = 0.0f;
        for(unsigned int i = 0; i < G_east.size(); i++){
            for(unsigned int j = 0; j < G_east.size(); j++){
                uH_ij += ipsi_east(j,i)*bvp.rbf.Value(pos_east[j],X,c2)*bvp.u.Value(pos_east[i],0.0f);
            }
        }
        return uH_ij;
    }
    std::cout << "Something went wrong computing H.";
    return 0.0f;
}

void Stencil::G_return(std::vector<int> & G_j, std::vector<double> & G){
	G_j.clear();
	G.clear();
	double norm = 1.0f/(counter_north + counter_south + counter_east + counter_west);
	if(pos_north.size() > 0){
        //norm = 1.0f/(counter_north);
		for(unsigned int i = 0; i < index_north.size(); i++){
			G_j.push_back(index_north[i]);
			G.push_back(G_north[i]*norm);
		}
	}
	if(pos_south.size() > 0){
        //norm = 1.0f/(counter_south);
		for(unsigned int i = 0; i < index_south.size(); i++){
			G_j.push_back(index_south[i]);
			G.push_back(G_south[i]*norm);
		}
	}
	if(pos_east.size() > 0){
        //norm = 1.0f/(counter_east);
		for(unsigned int i = 0; i < index_east.size(); i++){
			G_j.push_back(index_east[i]);
			G.push_back(G_east[i]*norm);
		}
	}
	if(pos_west.size() > 0){
        //norm = 1.0f/(counter_west);
		for(unsigned int i = 0; i < index_west.size(); i++){
			G_j.push_back(index_west[i]);
			G.push_back(G_west[i]*norm);
		}
	}

}
void Stencil::G_return_withrep(std::vector<int> & G_j, std::vector<double> & G, int N_tray){
    G_j.clear();
    G.clear();
    double norm = 1.0f/(N_tray*1.0f);
    std::map<int,double> G_map;
    if(pos_north.size() > 0){
        for(unsigned int i = 0; i < index_north.size(); i++){
            G_map[index_north[i]] = 0.0f;
        }
    }
    if(pos_south.size() > 0){
        for(unsigned int i = 0; i < index_south.size(); i++){
            G_map[index_south[i]] = 0.0f;
        }
    }
    if(pos_east.size() > 0){
        for(unsigned int i = 0; i < index_east.size(); i++){
            G_map[index_east[i]] = 0.0f;
        }
    }
    if(pos_west.size() > 0){
        for(unsigned int i = 0; i < index_west.size(); i++){
            G_map[index_west[i]] = 0.0f;
        }
    }
    if(pos_north.size() > 0){
        for(unsigned int i = 0; i < index_north.size(); i++){
            G_map[index_north[i]] =+ G_north[i]*norm;
        }
    }
    if(pos_south.size() > 0){
        for(unsigned int i = 0; i < index_south.size(); i++){
            G_map[index_south[i]] =+ G_south[i]*norm;
        }
    }
    if(pos_east.size() > 0){
        for(unsigned int i = 0; i < index_east.size(); i++){
            G_map[index_east[i]] =+ G_east[i]*norm;
        }
    }
    if(pos_west.size() > 0){
        for(unsigned int i = 0; i < index_west.size(); i++){
            G_map[index_west[i]] =+ G_west[i]*norm;
        }
    }
    for(std::map<int, double>::iterator it = G_map.begin();
        it != G_map.end();
        it ++){
        G_j.push_back(it->first);
        G.push_back(it->second);
    }

}
int Stencil::G_Test_update(Eigen::VectorXd X){
    if(AreSame(X(1),stencil_parameters[1])){
        //South stencil
        //std::cout << "\t south \n";
        for(unsigned int i = 0; i < G_south.size(); i++) G_south[i] = (double)index_south[i];
            return 1;
    }
    if(AreSame(X(1),stencil_parameters[3])){
        //North stencil
        for(unsigned int i = 0; i < G_north.size(); i++) G_north[i] = (double) index_north[i];
            return 1;
    }
    if(AreSame(X(0),stencil_parameters[0])){
        //West stencil
        //std::cout << "\t west \n";
        for(unsigned int i = 0; i < G_west.size(); i++) G_west[i] = (double) index_west[i];
            return 1;
    }
    if(AreSame(X(0),stencil_parameters[2])){
        //East stencil
        for(unsigned int i = 0; i < G_east.size(); i++) G_east[i] = (double) index_east[i];
            return 1;
    }
    std::cout << "Something went wrong computing H.";
    return 0;
}

void Stencil::G_Test_return(std::vector<int> & stencil_index, std::vector<double> & G){
    stencil_index.clear();
    G.clear();
    if(pos_north.size() > 0){
        //norm = 1.0f/(counter_north);
        for(unsigned int i = 0; i < index_north.size(); i++){
            stencil_index.push_back(index_north[i]);
            G.push_back(G_north[i]);
        }
    }
    if(pos_south.size() > 0){
        //norm = 1.0f/(counter_south);
        for(unsigned int i = 0; i < index_south.size(); i++){
            stencil_index.push_back(index_south[i]);
            G.push_back(G_south[i]);
        }
    }
    if(pos_east.size() > 0){
        //norm = 1.0f/(counter_east);
        for(unsigned int i = 0; i < index_east.size(); i++){
            stencil_index.push_back(index_east[i]);
            G.push_back(G_east[i]);
        }
    }
    if(pos_west.size() > 0){
        //norm = 1.0f/(counter_west);
        for(unsigned int i = 0; i < index_west.size(); i++){
            stencil_index.push_back(index_west[i]);
            G.push_back(G_west[i]);
        }
    }

}

void Stencil::Print(int node_index){
    char fname[100];
    FILE *pf;
    sprintf(fname,"Output/Debug/stencil_%d.txt", node_index);
    pf = fopen(fname, "w");
    fprintf(pf,"index,x,y,sten\n");
    for(unsigned int i = 0; i < index_south.size(); i++){
        fprintf(pf,"%d,%.3f,%.3f,south\n",index_south[i], pos_south[i][0], pos_south[i][1]);
    }
    for(unsigned int i = 0; i < index_east.size(); i++) fprintf(pf,"%d,%.3f,%.3f,east\n",index_east[i], pos_east[i](0), pos_east[i](1));
    for(unsigned int i = 0; i < index_north.size(); i++) fprintf(pf,"%d,%.3f,%.3f,north\n",index_north[i], pos_north[i](0), pos_north[i](1));
    for(unsigned int i = 0; i < index_west.size(); i++) fprintf(pf,"%d,%.3f,%.3f,west\n",index_west[i], pos_west[i](0), pos_west[i](1));   
    /*for(unsigned int i = 0; i < index_south.size(); i++) fprintf(pf,"%d\n",index_south[i]);
    for(unsigned int i = 0; i < index_east.size(); i++) fprintf(pf,"%d\n",index_east[i]);
    for(unsigned int i = 0; i < index_north.size(); i++) fprintf(pf,"%d\n",index_north[i]);
    for(unsigned int i = 0; i < index_west.size(); i++) fprintf(pf,"%d\n",index_west[i]);*/
    fclose(pf);
    sprintf(fname,"Output/Debug/boundary_stencil_%d.txt", node_index);
    pf = fopen(fname, "w");
    fprintf(pf,"x,y\n");
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[1]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[2], stencil_parameters[1]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[2], stencil_parameters[3]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[3]);
    fprintf(pf,"%.3f,%.3f\n",stencil_parameters[0], stencil_parameters[1]);
    fclose(pf);

}

