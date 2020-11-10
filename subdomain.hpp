#ifndef SUBDOMAIN
#define SUBDOMAIN
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>
#include "interface.hpp"
#include "BVP.hpp"
class Subdomain{
    private:
    //North, east, south and west interfaces
    Interface north, south, east, west;
    //Mesh of points inside the subdomain
    std::vector<Eigen::VectorXd> mesh;
    //Solution of the nodes of the Mesh
    std::vector<double> solution;
    public:
    //True if the subdomain is solved
    bool solved;
    //Label of the subdomain
    std::vector<int> label;
    //Default initialization
    Subdomain(void);
    //Initialization of the variables of an instance
    void Init(std::vector<int> subdomain_label, Interface north_interface, 
        Interface south_interface, Interface east_interface, Interface west_interface);
    //Solves the BVP in the subdomain
    void Solve(BVP bvp);
    //Prints the solution in a given filename
    void Print_Solution(char* filename);

};
#endif