#include "subdomain.hpp"
Subdomain::Subdomain(void){
    solved = false;
}
void Subdomain::Init(std::vector<int> subdomain_label, Interface north_interface, 
    Interface south_interface, Interface east_interface, Interface west_interface){
    label = subdomain_label;
    north = north_interface;
    south = south_interface;
    east = east_interface;
    west = west_interface;
}
//Solves the BVP in the subdomain
void Solve(BVP bvp){
    printf("To be coded\n");
}
//Prints the solution in a given filename
void Print_Solution(char* filename){
    printf("To be coded\n");
}