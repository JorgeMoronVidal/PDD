#ifndef STENCIL 
#define STENCIL
#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <fstream>
#include <algorithm>
#include "boundary.hpp"
#include "BVP.hpp"
#define EPSILON 1e-6
enum direction {North, South, East, West};
class Stencil
{
    private:
        //Inverse psi matrix for all the four stencils
	    Eigen::MatrixXd ipsi_north, ipsi_south, ipsi_east, ipsi_west;
        //Position of the nodes of the stencils
	    std::vector<Eigen::VectorXd> pos_north, pos_south, pos_east, pos_west;
        //2D indexes of the interfaces that originate the stencils 
        std::vector<int> index_north, index_south, index_east, index_west;
        //Cumulatives of the G matrix values associated to the different stencils
        std::vector<double> G_north, G_south, G_east, G_west;
        //Counter of the trayectoies that ends in each stencil
        uint32_t counter_north, counter_south, counter_east, counter_west;
        //Computes the inverse matrix of the psi matrix
        Eigen::MatrixXd Compute_ipsi(std::vector<Eigen::VectorXd> & sten_position, BVP bvp, double c2);
       
    public:
        //Parameters of the stencil boundary
        double stencil_parameters[4];
        //Parameters of the boundary of the problem
        double global_parameters[4];
        //Deafault initializer
	    Stencil(void);
        /*double comparison with epsilon tolerance
        double a One of the double numbers to compare
        double b the ther double number
        */
        bool AreSame(double a, double b);
        /*Proper intializer of the Stencil class when the node belongs to 4 interfaces 
        std::vector<int> interface_east East interface the node belongs to
        std::vector<int> interface_south South interface the node belongs to
        std::vector<int> interface_west West interface the node belongs to
        std::vector<int> interface_north North interface the node belongs to 
        std::vector<int> current_index Index of the interface
        std::vector<Eigen::VectorXd> & node_position Position of the nodes in the problem
        std::vector<std::vector<int> > & node_index Indexes of the nodes in the problem
        std::map<std::vector<int>, int> & interface_map Map which maps the interfaces 2D index into the interfaces 1D index
        std::vector<int> & number_interfaces Number of interfaces in each direction
        double g_parameters Parameters of the stencil boudnary 
        */
        void Init(std::map<direction, std::vector<int> > s_index,
                  std::map<direction, std::vector<double> > s_x,
                  std::map<direction, std::vector<double> > s_y,
                  double *g_parameters);
        /*Resets the stencil object*/
        void Reset(void);
        /*Computes the inverse psi matrices given
        BVP bvp The boundary value problem 
        double c2 The constan for the mesheless interpolator
        */
        void Compute_ipsi(BVP bvp, double c2);
        /* Updates G matrix 
        Eigen::VectorXd X position where the trayectory ended
        double Y Y FKAC value where the trayectoy ended
        BVP bvp Boundary value problem
        double c2 Constant for the mesheless interpolator*/
        int G_update(Eigen::VectorXd X, double Y, BVP bvp, double c2);
        /* Updates G matrix for the test fucntion
        Eigen::VectorXd X position where the trayectory ended*/
        int G_Test_update(Eigen::VectorXd X);
        /* Test the Meshless interpolator in the stencils
        Eigen::VectorXd X position where the trayectory ended
        BVP bvp Boundary value problem
        double c2 Constant for the mesheless interpolator*/
        double Test_Interpolator(Eigen::VectorXd X, BVP bvp, double c2);
        /* Return G matrix if no nodes are repeated in the stencil
        std::vector<int> & stencil_index The indexes of the nodes that compound the stencil
        std::vector<double> & G The values of the G matrix
        */
        void G_return(std::vector<int> & stencil_index, std::vector<double> & G);
        /* Return G matrix taking care of the repeated nodes of the stencil
        std::vector<int> & stencil_index The indexes of the nodes that compound the stencil
        std::vector<double> & G The values of the G matrix
        int N_tray Total number of trayectories 
        */
        void G_return_withrep(std::vector<int> & stencil_index, std::vector<double> & G, int N_tray);
        /*Returns G function when the meshless interpolator is being tested:
        std::vector<int> & stencil_index Indexes of the stencil nodes
        std::vector<double> & G Values od the G matrix
        */
        void G_Test_return(std::vector<int> & stencil_index, std::vector<double> & G);
        /*Prints the stencil parameters
        int node_index index of the node whose stencil we are trying to compute
        */
        void Print(int node_index);
};
#endif