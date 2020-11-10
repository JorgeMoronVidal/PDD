#ifndef INTERFACE
#define INTERFACE
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>
struct Interface{
    //True if all the nodes of the interface are solved 
    bool solved; 
    //Check if the solution of a given node has been computed
    std::vector<bool> solved_nodes;
    //2D label of the interface
    std::vector<int> label;
    //Indexes of the nodes of the interface
    std::vector<int> index;
    //Positions of the nodes of the interface
    std::vector<Eigen::VectorXd> position;
    //Solution of the nodes of the interface
    std::vector<double> solution;
    //Initialization of the variables of an instance
    void Init(std::vector<int> lab, std::vector<int> ind);
    //Returns true if the node with index ind is in the interface
    bool In (int ind);
    //North interface or interfaces
    std::vector<std::vector<int> > N(void);
    //South interface or interfaces
    std::vector<std::vector<int> > S(void);
    //East interface or interfaces
    std::vector<std::vector<int> > E(void);
    //West interface or interfaces
    std::vector<std::vector<int> > W(void);
    //Returns the position of a node 
    Eigen::VectorXd Node_Position(int node_index);
    //Print Interface 
    void Print(int interface_i);
};
#endif