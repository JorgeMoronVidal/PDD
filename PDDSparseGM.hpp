//mpiexec --mca orte_base_help_aggregate 0  -np 4 ./main
#ifndef PDDSPARSEGM
#define PDDSPARSEGM
//Fraction of the interface the stencil is elongated
#define STEN_ELONG 0.5
//MPI TAGS
#define REQUEST_NODE 1
#define REPLY_NODE 2
#define NODE_POSITION 3
#define PARAMETERS 4
#define SIZES 5
#define IND_NORTH 6
#define IND_SOUTH 7
#define IND_EAST 8
#define IND_WEST 9
#define X_NORTH 10
#define X_SOUTH 12
#define X_EAST 14
#define X_WEST 16
#define Y_NORTH 11
#define Y_SOUTH 13
#define Y_EAST 15
#define Y_WEST 17
#define TAG_Gi 20
#define TAG_Gj 21
#define TAG_Gval 22
#define TAG_Bi 30
#define TAG_Bval 31
#include "GMSolver.hpp"
#include "interface.hpp"
#include "stencil.hpp"
#include "subdomain.hpp"
#include <mpi.h>
#include <fstream>

typedef Eigen::Triplet<double,int> T;

class PDDSparseGM{
    friend class GMSolver;
    private:
        /*
        MPI-related variables
        */
        MPI_Comm world, workers;
        MPI_Group world_group, worker_group;
        MPI_Status status;
        int ranks[1],work_control[2],numprocs, myid, server, workerid;
        /*
        Interfaces of the problem
        */
        std::vector<Interface> interfaces;
        /*
        Subdomains of the problem
        */
        std::vector<Subdomain> subdomains;
        /*-maps each interface 2D index in its 1D index*/
        std::map<std::vector<int>, int> interface_map;
        /*
         -Number of interfaces (And subdomains) per direction.
         -Nuber of nodes per interface.
        */
        std::vector<int> iN,nN;
        /*File were the parameters of the problem are written*/
        std::string configuration_file;
        /*-h0 is the value of time discretization for the lowest level.
          -T_start is the time value where we want to compute dhe solution
          -eps represents the rse objective of the stimator
          -parameters stores the boundary parameters
          -factor that multiplies the internodal distance
          -constant of the RBF interpolator 
          */
        double h0, T_start, eps, fac, c2;
        std::vector<double> parameters;
        /*-SW is the south west point of the domain.
          -NE is the north east point of the domain.
        */
       Eigen::VectorXd SW, NE;
       /*G and B storage vector*/
       std::vector<double> G, B;
       std::vector<int>  G_j, G_i, B_i;
       /*Triplet's vector*/
       std::vector<T> T_vec_G, T_vec_B;
        /*
          -N is the initial number of trayectories 
        */
       int N;
       /*Stencil of the node*/
       Stencil node_stencil;
       /*Position of the Node being solver*/
       Eigen::VectorXd position;
       /*Starts the MPI values and structures*/
       void MPI_Configuration(int argc, char *argv[]);
       /*Interface, subdomain and node structure intialization*/
       void Fullfill_interfaces(void);
       /*Returns the interfaces the node belongst to*/
       std::vector<std::vector<int> > Get_Interfaces(int index);
       /*Returns the stencils labels*/
       std::map<direction, std::vector<std::vector <int> > > Labels_Stencil(int index);
       /*Send node information*/
       void Send_Node(int index);
       /*Receive node information*/
       bool Receive_Node(void);
       /*Sends the stencil data to the worker*/
       void Send_Stencil_Data(int index);
       /*Receives the stencil data from the server*/
       Stencil Recieve_Stencil_Data(void);
       /*Send G matrix and B Matrix*/
       void Send_G_B(void);
       /*Receive G matrix and B Matrix*/
       void Receive_G_B(void);
       /*Position of a given node*/
       Eigen::VectorXd Node_Position(int node_index);
       /*Prints metadata file*/
       void Print_metadata(void);
       /*Prints solution file*/
       void Print_solution(BVP bvp);
       /*Prints Domains files*/
       void Print_subdomain(BVP bvp);
    public:
       /*Reads the problem information from a file.*/
       void ReadFromFile(std::string file);
       /*Initialization of the class*/
       PDDSparseGM(int argc, char *argv[]);
       /*Initialization of the class with a configuration file*/
       PDDSparseGM(int argc, char *argv[], std::string file);
       /*Solves the problem given by the configuration file*/
       void Solve(BVP bvp, std::string file);
       /*Solves the problem given by the configuration file*/
       void Solve(BVP bvp);
       /*Prints the problem parameters on screen*/
       void Print_Problem(void); 
       /*Prints the content of the Interface vector on screen*/
       void Print_Interface(void);
};
#endif