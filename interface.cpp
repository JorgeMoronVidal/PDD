#include "interface.hpp"
void Interface::Init(std::vector<int> lab, std::vector<int> ind){
    label = lab;
    ind = index;
}

bool Interface::In(int node_index){
    for(std::vector<int>::iterator ind = index.begin();
    ind != index.end();
    ind ++){
        if(*ind == node_index) return true;
    }
    return false;
}

std::vector<std::vector<int> > Interface::N(void){
    std::vector<int> aux;
    std::vector<std::vector<int> > output;
    if(label[1]%2 == 0){
        //Vertical Interface
        aux = label;
        aux[1] ++;
        output.push_back(aux);
        aux[0] ++;
        output.push_back(aux);

    } else {
        //Horizontal Interface
        aux = label;
        aux[1] += 2;
        output.push_back(aux);
    }
    return output;
}
std::vector<std::vector<int> > Interface::S(void){
    std::vector<int> aux;
    std::vector<std::vector<int> > output;
    if(label[1]%2 == 0){
        //Vertical Interface
        aux = label;
        aux[1] -= 1;
        output.push_back(aux);
        aux[0] ++;
        output.push_back(aux);

    } else {
        //Horizontal Interface
        aux = label;
        aux[1] -= 2;
        output.push_back(aux);
    }
    return output;
}
std::vector<std::vector<int> > Interface::E(void){
    std::vector<int> aux;
    std::vector<std::vector<int> > output;
    if(label[1]%2 == 0){
        //Vertical Interface
        aux = label;
        aux[0] ++;
        output.push_back(aux);

    } else {
        //Horizontal Interface
        aux = label;
        aux[1] --;
        output.push_back(aux);
        aux[1] += 2;
        output.push_back(aux);
    }
    return output;
}
std::vector<std::vector<int> > Interface::W(void){

    std::vector<int> aux;
    std::vector<std::vector<int> > output;
    if(label[1]%2 == 0){
        //Vertical Interface
        aux = label;
        aux[0] --;
        output.push_back(aux);

    } else {
        //Horizontal Interface
        aux = label;
        aux[0] --;
        aux[1] --;
        output.push_back(aux);
        aux[1] += 2;
        output.push_back(aux);
    }
    return output;
}
Eigen::VectorXd Interface::Node_Position(int node_index){
    Eigen::VectorXd output;
    output.resize(2);
    for(unsigned int i = 0; i < index.size(); i++){
        if(index[i] == node_index) output = position[i];
    }
    return output;
}
void Interface::Print(int interface_i){
    char filename[100];
    FILE *pFile;
    sprintf(filename,"Output/Debug/Interface_%d.txt", interface_i);
    pFile = fopen(filename, "w");
    fprintf(pFile,"index,x,y,solved,solution\n");
    for(int i = 0; i < (int) index.size(); i++){
        fprintf(pFile,"%d\t %f \t %f \t %d \t %f \n",
        index[i], position[i][0], position[i][1],
        solved_nodes[i], solution[i]);
    }
    fclose(pFile);
}