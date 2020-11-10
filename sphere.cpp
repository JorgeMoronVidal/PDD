#include "sphere.hpp"
#include "iostream"
double Sphere(double* params, 
            Eigen::VectorXd & position, 
            Eigen::VectorXd & exitpoint,
            Eigen::VectorXd & normal){
    
    //Give the distance to a sphere of radious parameters[0] centered in 
    //(parameters[1],...,parameters[N])
    double r = 0.0f; 

    for(int i = 0; i < position.size(); i++){
        
        r += pow(position(i)-params[i+1],2.0);

    }

    r = sqrt(r);
    normal = position/r;
    for(int i = 0; i < position.size(); i++){
        exitpoint(i) = normal(i) * params[0];
    }
    
    if(r < params[0]){
        return r - params[0];
    }

    return r - params[0];
    
}
            
void Plot_Sphere(double * params, Eigen::VectorXd & position){
    FILE *f;
    double N = 200.0;
    if(position.size() == 2){
        f = fopen("Output/trayectories/surface.txt", "w");
        fprintf(f,"X,Y\n");
        for(double phi = 0 ; phi <= 2 * M_PI+0.001 ; phi += (2.0/N) * M_PI){
            fprintf(f,"%e,%e\n",params[1] + params[0]*cos(phi), params[2] + params[0]*sin(phi));
        }
        fclose(f);
    } else {
        if (position.size() == 3){
            f = fopen("Output/trayectories/surface.txt", "w");
            fprintf(f,"X,Y,Z\n");
            for(double theta = 0 * M_PI; theta <= M_PI+0.01; theta += (1.0f/N) * M_PI){
                for(double phi = 0 * M_PI; phi <= 2 * M_PI+0.001; phi += (1.0f/(0.5*N)) * M_PI){
                    fprintf(f,"%e,%e,%e\n",params[1] + params[0]*sin(theta)*cos(phi), params[2] + params[0]*sin(theta)*sin(phi), params[3] + params[0]*cos(theta));
                }
            }
            fclose(f);
        } else {
            printf("Is not possible to plot surfaces with more than 3 dimensions. \n");
        }
    }
}