#include "PDDSparseGM.hpp"
//mpiexec -np 4 xterm -e gdb ./main -ex run
int main(int argc, char *argv[]){
    std::string config("configuration.txt");
    PDDSparseGM PDDS(argc,argv,config);
    BVP bvp;
    PDDS.Solve(bvp);
    return 0;
}