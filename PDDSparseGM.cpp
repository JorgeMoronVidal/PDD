#include "PDDSparseGM.hpp"
void PDDSparseGM::MPI_Configuration(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &numprocs);
    MPI_Comm_rank(world, &myid);
    server = numprocs-1;
    MPI_Comm_group(world, &world_group);
    ranks[0] = server;
    //We exclude a group in world to create workers
    MPI_Group_excl(world_group, 1, ranks, &worker_group);
    //Workers comunicator is created
    MPI_Comm_create(world, worker_group, &workers);
    //We free both comunicators groups
    MPI_Group_free(&worker_group);
    MPI_Group_free(&world_group);
}
PDDSparseGM::PDDSparseGM(int argc, char *argv[]){
  iN.resize(2);
  nN.resize(2);
  h0 = 0.0;
  T_start =INFINITY;
  eps = 0.5;
  SW.resize(2);
  NE.resize(2);
  N = 0;
  fac = 1;
  MPI_Configuration(argc, argv);
}
PDDSparseGM::PDDSparseGM(int argc, char *argv[],std::string file){ 
  MPI_Configuration(argc, argv);
  ReadFromFile(file);
  parameters.resize(4);
  parameters[0] = SW[0];
  parameters[1] = SW[1];
  parameters[2] = NE[0];
  parameters[3] = NE[1];
  if(myid == server){
    Fullfill_interfaces();
  }
}
void PDDSparseGM::ReadFromFile(std::string file){
  std::ifstream myfile(file.c_str());
    std::string line, cent;
    //If T is not specified then the problem is suposed to be elliptic
    T_start = INFINITY;
    if (myfile)  // same as: if (myfile.good())
    {   while (getline( myfile, line )){  // same as: while (getline( myfile, line ).good())

            if (line.c_str()[0]!='%') //First of all we check if the line is a comment
            { 
                cent.clear(); //Centinel needs to be empty

                for(std::string::iterator it = line.begin(); it != line.end(); it++){
                //Loop over all the characters of the line
                    
                    if(*it != ' ' && *it != '\t'){
                    //If it is not a space nor a tabulation then we concatenate the term to cent
                        cent += *it;
                    }
                    else{
                    //If we get a space or a tabulation then we check the string
                        if(cent == "h0="){
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            h0 = atof(cent.c_str());
                            break;
                        }
                        if(cent == "eps="){
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            eps = atof(cent.c_str());
                            break;
                        }
                        if(cent == "fac="){
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            fac = atof(cent.c_str());
                            break;
                        }
                        if(cent == "T="){
                         
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            T_start = atof(cent.c_str());
                            break;
                        }
                        if(cent == "SW="){
                        //South West point of the domain 
                            SW.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                                while((*it != ' ') && (*it != '\t') && (it != line.end())){

                                cent += *it;
                                it++;

                                }
                                SW(i) = atof(cent.c_str());
                                cent.clear();
                                while( *it == ' ' or *it == '\t'){
                                        it++;
                                }
                    
                            }

                            break;
                        }

                        if(cent == "NE="){
                        //South West point of the domain 
                            NE.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                            while((*it != ' ') && (*it != '\t') && (it != line.end())){

                                cent += *it;
                                it++;

                            }
                            NE(i) = atof(cent.c_str());
                            cent.clear();
                            while( *it == ' ' or *it == '\t'){
                                        it++;
                            }
                    
                            }

                            break;
                        }

                        if(cent == "N="){
                        //Number of trayectories per node

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            N = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_subdo="){
                        //Number of domains in each direction

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            iN.resize(2);
                            iN[0] = atoi(cent.c_str());
                            iN[1] = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_node="){
                        //Number nodes per interface

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            nN.resize(2);
                            nN[0] = atoi(cent.c_str());
                            nN[1] = atoi(cent.c_str());
                            break;
                        }
                        
                    }
                }
            }
        }
        myfile.close();

    }

    else {
    
        std::cout << "WARNING: Cofiguration file was not found\n";
        std::cout << "Make sure it exist a configuration.txt file in the program's directory.\n";
    }    
}
void PDDSparseGM::Fullfill_interfaces(void){
   /*Initial values, long and  and shot increments for each direction*/
   std::vector<Eigen::VectorXd> P_cent, lincrement, sincrement, start, end;
   Eigen::VectorXd vaux;
   std::vector<std::vector<int>> node_index;
   interfaces.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   subdomains.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   node_index.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   vaux.resize(SW.size());

   //Increments are defined
   //ShortIncrement
   for (int i = 0; i < SW.size() ; i++){
        for(int j = 0; j < vaux.size(); j++){
            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0* ((nN[j])*iN[j]));

            } else {

                vaux(j) = 0.0;

            }
        }
        sincrement.push_back(vaux);
    }
    //LongIncrement
    for (int i = 0; i < SW.size() ; i++){

        for(int j = 0; j < vaux.size(); j++){

            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0* (iN[j]));

            } else {

                vaux(j) = 0.0;

            }
        }

        lincrement.push_back(vaux);
    }
    //Initial Point
   
    for (int i = 0; i < SW.size(); i++){
            vaux = SW;
            P_cent.push_back(vaux);
    }
    //c2 is computed as the square of fac*internodal space.
    c2 = pow(fac*sincrement[0].norm(),2.0);
    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];
    P_cent[1] += lincrement[0];

    int cen = 0;
    for (int i = 0; i < iN[1] - 1; i++){

        for (int j = 0; j < iN[0]; j++){

            start.push_back(P_cent[0] + lincrement[0]*j);
            end.push_back(P_cent[0] + lincrement[0]*(j+1));

            interfaces[cen].label.push_back(j);
            interfaces[cen].label.push_back(2*i + 1);
            subdomains[cen].label.push_back(j);
            subdomains[cen].label.push_back(2*i + 1);
            cen++;

            start.push_back(P_cent[1] + lincrement[1]*j);
            end.push_back(P_cent[1] + lincrement[1]*(j+1));

            interfaces[cen].label.push_back(i);
            interfaces[cen].label.push_back(2*j);
            subdomains[cen].label.push_back(i);
            subdomains[cen].label.push_back(2*j);

            cen++;
        }


        P_cent[0] += lincrement[1];
        P_cent[1] += lincrement[0];

    }

    for(int i = 0; i < (int)interfaces.size(); i++){
        interface_map[interfaces[i].label] = i;
    }
    
    //stvaux stores the current interface index.
    //cross stores the index of the nodes that constitutes 
    std::vector<int> stvaux, cross;
    unsigned int index, cent = 0, cross_cent = 0;
    cross.resize((iN[0]-1)*(iN[1]-1));
    //As It is done in Paco's original code we start by fullfilling the vertical interfaces
    cross.clear();
    for(int i = 0; i < iN[0] - 1; i++){
        for(int j = 0; j <= iN[1]*2 - 2; j += 2){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            node_index[index].push_back(cent);
            if(j != 0) cross[(j/2-1)*(iN[1]-1) + i] = cent;
            for(int k = 0; k < nN[1]; k ++){
                cent ++;
                node_index[index].push_back(cent);
            }
            if(j == iN[1]*2 -2) {
                cent ++;
            }
        }
    }
    //Horizontal interfaces
    for(int j = 1; j <= 2*iN[1]-3; j +=2 ){
        for(int i = 0; i < iN[0]; i ++){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            if(i != 0){
                node_index[index].push_back(cross[cross_cent]);
                cross_cent ++;
            } else {
                node_index[index].push_back(cent);
                cent ++;
            }
            for(int k = 0; k < nN[0] - 1; k ++){
                node_index[index].push_back(cent);
                cent ++;
            }
            if(i != iN[0] -1 ) {
                node_index[index].push_back(cross[cross_cent]);
            } else {
                node_index[index].push_back(cent);
                cent ++;
            }
        }
    }
    nNodes = (int)cent;
    std::cout << "nNodes is " << nNodes << std::endl;
    //Loop over the points in the system
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        index = interface_map[(*it).label];
        (*it).index.resize(node_index[index].size());
        (*it).solution.resize(node_index[index].size());
        (*it).position.resize(node_index[index].size());
        sincrement[0] = (end[index]-start[index])/(node_index[index].size()-1);
        vaux = start[index];
        for(int i = 0; i < (int) node_index[index].size(); i++){
            (*it).index[i] = node_index[index][i];
            (*it).position[i] = vaux;
            vaux += sincrement[0];
        }
    }
}
void PDDSparseGM::Solve(BVP bvp){
    bool done=true;
    if(myid==server) Print_Problem();
    //Compute the PDDSparse Matrix
    if(myid==server){
        FILE * pFile;
        pFile = fopen ("Output/Debug/Node_position.txt","w");
        fprintf(pFile, "index,x,y\n");
        fclose(pFile);
        //Server distributes the work between the workers process
        for(int node_index = 0; node_index < nNodes; node_index++){
            work_control[1] = 1;
            Send_Node(node_index);
            Send_Stencil_Data(work_control[0]);
        }
        //G and B are received from the workers
        for(int process = 0; process < server; process++){
             Receive_G_B();
        }
        Compute_Solution(bvp);
        printf("Process %d ended its work.\n",myid);
    } else {
        double start = MPI_Wtime();
        //G and B storage vectors for each node
        double B_temp;
        std::vector<double> G_temp;
        std::vector<int>  G_j_temp;
        //Stencil
        Stencil stencil;
        GMSolver solver(bvp ,parameters, h0, (unsigned int) myid +1);
        do{
            done = Receive_Node();
            if(!done){
                G_temp.clear();
                B_temp = 0.0;
                G_j_temp.clear();
                G_temp.clear();
                stencil = Recieve_Stencil_Data();
                stencil.Compute_ipsi(bvp,c2);
                solver.Solve(position,c2,stencil,G_j_temp,G_temp,B_temp,N);
                B.push_back(B_temp);
                B_i.push_back(work_control[0]);
                for(unsigned int k = 0;k < G_temp.size(); k ++){
                    G_i.push_back(work_control[0]);
                    G_j.push_back(G_j_temp[k]);
                    G.push_back(G_temp[k]);
                }
            }
        }while(!done);
        printf("Process %d is done solving nodes \n", myid);
        Send_G_B();
        printf("Process %d is done sending G \n", myid);
        double end = MPI_Wtime();
        printf("Process %d used %f  hours \n", myid, (start-end)/3600);
    }
    MPI_Finalize();
}
void PDDSparseGM::Solve(BVP bvp, std::string file){

    ReadFromFile(file);
    if(myid == server){
          Fullfill_interfaces();
    }
    Solve(bvp);
}
void PDDSparseGM::Send_Node(int index){
    int work_control_aux[2];
    MPI_Recv(work_control_aux, 2, MPI_INT, MPI_ANY_SOURCE, REQUEST_NODE, world, &status);
    work_control[0] = index;
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    double pos[2];
    pos[0] = Node_Position(index)[0];
    pos[1] = Node_Position(index)[1];
    MPI_Send(pos, 2, MPI_DOUBLE, status.MPI_SOURCE, NODE_POSITION, world);
    if(work_control[1] == 1){
        FILE *pFile;
        pFile = fopen ("Output/Debug/Node_position.txt","a");
        fprintf(pFile, "%d,%.4f,%.4f \n",index,pos[0],pos[1]);
        fclose(pFile);
        printf("Sending node %d.\n",work_control[0]);
        Send_Stencil_Data(index);
    } else {
        printf("Ending process.\n");
    }
}
std::vector<std::vector<int> > PDDSparseGM::Get_Interfaces(int index){
    std::vector<std::vector<int> > output;
    //printf("Node %d is in interfaces: ",index);
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        if((*it).In(index) == true) {
            output.push_back((*it).label);
            //printf("[%d %d]", (*it).label[0],(*it).label[1]);
        }
    }
    //printf("\n");
    return output;
}
std::map<direction, std::vector<std::vector <int> > > PDDSparseGM::Labels_Stencil(int index){
    std::map<direction, std::vector<std::vector <int> > > output;
    std::vector< std::vector<int> > node_interfaces;
    Interface inter;
    std::vector< int > interface_east, interface_west,interface_north, interface_south;
    int max_v = 0, max_h = 0;
    node_interfaces = Get_Interfaces(index);
    switch(node_interfaces.size()){
        case 1:
            inter = interfaces[interface_map[node_interfaces[0]]];
            for(int i = 0; i < (int) inter.N().size(); i ++){
                if(interface_map.count(inter.N()[i]) > 0) output[North].push_back(inter.N()[i]);
            }
            for(int i = 0; i < (int) inter.S().size(); i ++){
                if(interface_map.count(inter.S()[i]) > 0) output[South].push_back(inter.S()[i]);
            }
            for(int i = 0; i < (int) inter.E().size(); i ++){
                if(interface_map.count(inter.E()[i]) > 0) output[East].push_back(inter.E()[i]);
            }
            for(int i = 0; i < (int) inter.W().size(); i ++){
                if(interface_map.count(inter.W()[i]) > 0) output[West].push_back(inter.W()[i]);
            }
            break;
        case 4:
            max_h = std::max(node_interfaces[0][0], node_interfaces[1][0]);
            max_h = std::max(max_h, node_interfaces[2][0]);
            max_h = std::max(max_h, node_interfaces[3][0]);
            max_v = std::max(node_interfaces[0][1], node_interfaces[1][1]);
            max_v = std::max(max_v, node_interfaces[2][1]);
            max_v = std::max(max_v, node_interfaces[3][1]);
            if(node_interfaces[0][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[0][0] == max_h) interface_east = node_interfaces[0];
                else interface_west = node_interfaces[0];
            } else {
                //Vertical interface
                if(node_interfaces[0][1] == max_v) interface_north = node_interfaces[0];
                else interface_south = node_interfaces[0];
            }
            if(node_interfaces[1][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[1][0] == max_h) interface_east = node_interfaces[1];
                else interface_west = node_interfaces[1];
            } else {
                //Vertical interface
                if(node_interfaces[1][1] == max_v) interface_north = node_interfaces[1];
                else interface_south = node_interfaces[1];
            }
            if(node_interfaces[2][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[2][0] == max_h) interface_east = node_interfaces[2];
                else interface_west = node_interfaces[2];
            } else {
                //Vertical interface
                if(node_interfaces[2][1] == max_v) interface_north = node_interfaces[2];
                else interface_south = node_interfaces[2];
            }
            if(node_interfaces[3][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[3][0] == max_h) interface_east = node_interfaces[3];
                else interface_west = node_interfaces[3];
            } else {
                //Vertical interface
                if(node_interfaces[3][1] == max_v) interface_north = node_interfaces[3];
                else interface_south = node_interfaces[3];
            }
            inter = interfaces[interface_map[interface_west]];
            if(interface_map.count(inter.N()[0]) > 0) output[North].push_back(inter.N()[0]);
            inter = interfaces[interface_map[interface_east]];
            if(interface_map.count(inter.N()[0]) > 0) output[North].push_back(inter.N()[0]);
            inter = interfaces[interface_map[interface_west]];
            if(interface_map.count(inter.S()[0]) > 0) output[South].push_back(inter.S()[0]);
            inter = interfaces[interface_map[interface_east]];
            if(interface_map.count(inter.S()[0]) > 0) output[South].push_back(inter.S()[0]);
            inter = interfaces[interface_map[interface_south]];
            if(interface_map.count(inter.E()[0]) > 0) output[East].push_back(inter.E()[0]);
            inter = interfaces[interface_map[interface_north]];
            if(interface_map.count(inter.E()[0]) > 0) output[East].push_back(inter.E()[0]);
            inter = interfaces[interface_map[interface_south]];
            if(interface_map.count(inter.W()[0]) > 0) output[West].push_back(inter.W()[0]);
            inter = interfaces[interface_map[interface_north]];
            if(interface_map.count(inter.W()[0]) > 0) output[West].push_back(inter.W()[0]);
            break;
        default: 
            printf("Something went wrong in Labels_Stencil.\n");
    }
    /*
    printf("South: ");
    for(int i = 0; i < (int) output[South].size(); i++) printf("[%d %d] ",output[South][i][0],output[South][i][1]);
    printf("\n");
    printf("North: ");
    for(int i = 0; i < (int) output[North].size(); i++) printf("[%d %d] ",output[North][i][0],output[North][i][1]);
    printf("\n");
    printf("East: ");
    for(int i = 0; i < (int) output[East].size(); i++) printf("[%d %d] ",output[East][i][0],output[East][i][1]);
    printf("\n");
    printf("West: ");
    for(int i = 0; i < (int) output[West].size(); i++) printf("[%d %d] ",output[West][i][0],output[West][i][1]);
    printf("\n"); */
    return output;
}
void PDDSparseGM::Send_Stencil_Data(int index){
    std::map<direction, std::vector<std::vector <int> > > labels = Labels_Stencil(index);
    //The parameters of the stencil are computed
    double s_params[4];
    for(int i = 0; i < 4; i++) s_params[i] = parameters[i];
    if(labels[South].size()>0){
        s_params[0] = interfaces[interface_map[labels[South][0]]].position[0][0];
        s_params[1] = interfaces[interface_map[labels[South][0]]].position[0][1];
    } else {
        if(labels[West].size()>0){
        s_params[0] = interfaces[interface_map[labels[West][0]]].position[0][0];
        s_params[1] = interfaces[interface_map[labels[West][0]]].position[0][1];
        }
    }
    int aux_size;
    if(labels[North].size()>0){
        aux_size = (int) interfaces[interface_map[labels[North][labels[North].size()-1]]].position.size();
        s_params[2] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][0];
        s_params[3] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][1];
    } else {
        if(labels[East].size()>0){
        aux_size = (int) interfaces[interface_map[labels[East][labels[East].size()-1]]].position.size();
        s_params[2] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][0];
        s_params[3] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][1];
        }
    }
    //The indexes and positions of the stencils are computed
    std::map<direction, std::vector<int> > s_index;
    std::map<direction, std::vector<double> > s_x, s_y;
    std::vector<direction> s_dir;
    Interface s_inter, aux_inter;
    std::vector<int> aux_label;
    s_dir.resize(4);
    s_dir[0] = North; s_dir[1] = South; s_dir[2] = East; s_dir[3] = West;
    for(std::vector<direction>::iterator it = s_dir.begin(); it != s_dir.end(); it++){
    if(labels[*it].size() > 0){
        s_inter = interfaces[interface_map[labels[*it][0]]];
        aux_label.resize(2);
        aux_label[0] = s_inter.label[0]; aux_label[1] = s_inter.label[1];
        switch(*it){
            case North:
                aux_label[0] --;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case South:
                aux_label[0] --;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
                break;
            case East:
                aux_label[1] += -2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case West:
                aux_label[1] += -2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            default:
                printf("Something went wrong during stencil formation. 1\n");
        }
        for(int i = 0; i < (int) labels[*it].size(); i++){
            s_inter = interfaces[interface_map[labels[*it][i]]];
            for(int j = 0;
            j < (int) s_inter.index.size(); j ++){
                s_index[*it].push_back(s_inter.index[j]);
                s_x[*it].push_back(s_inter.position[j][0]);
                s_y[*it].push_back(s_inter.position[j][1]);
            }
        }
        aux_label.resize(2);
        aux_label[0] = s_inter.label[0]; aux_label[1] = s_inter.label[1];
        switch(*it){
            case North:
                aux_label[0] ++;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case South:
                aux_label[0] ++;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
                break;
            case East:
                aux_label[1] += 2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                std::cout << "East nodes " << (int) ((STEN_ELONG)*(aux_inter.index.size())) << std::endl;
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case West:
                aux_label[1] += 2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            default:
                printf("Something went wrong during stencil formation. 2\n");
        }
    }
    }
    //Parameters of the stencils are sent to the worker
    MPI_Send(s_params, 4, MPI_DOUBLE, status.MPI_SOURCE, PARAMETERS, world);
    int sizes[4];
    sizes[0] = (int) s_index[North].size();sizes[1] = (int) s_index[South].size();
    sizes[2] = (int) s_index[East].size();sizes[3] = (int) s_index[West].size();
    //Sizes of the stencil's vectors are sent to the worker
    MPI_Send(sizes, 4, MPI_INT, status.MPI_SOURCE, SIZES, world);  
    //Stencil indexes and positions are sent to the worker 
    int *i_north, *i_south, *i_west, *i_east;
    double *x_north, *x_south, *x_west, *x_east, 
           *y_north, *y_south, *y_west, *y_east;
    //North stencil
    i_north = new int[sizes[0]];
    x_north = new double[sizes[0]];
    y_north = new double[sizes[0]];
    for(int i = 0; i < sizes[0]; i++){
        i_north[i] = s_index[North][i];
        x_north[i] = s_x[North][i];
        y_north[i] = s_y[North][i];
    }
    MPI_Send(i_north, sizes[0], MPI_INT, status.MPI_SOURCE, IND_NORTH, world);
    MPI_Send(x_north, sizes[0], MPI_DOUBLE, status.MPI_SOURCE, X_NORTH, world);
    MPI_Send(y_north, sizes[0], MPI_DOUBLE, status.MPI_SOURCE, Y_NORTH, world);
    delete i_north; delete x_north; delete y_north;
    //South stencil
    i_south = new int[sizes[1]];
    x_south = new double[sizes[1]];
    y_south = new double[sizes[1]];
    for(int i = 0; i < sizes[1]; i++){
        i_south[i] = s_index[South][i];
        x_south[i] = s_x[South][i];
        y_south[i] = s_y[South][i];
    }
    MPI_Send(i_south, sizes[1], MPI_INT, status.MPI_SOURCE, IND_SOUTH, world);
    MPI_Send(x_south, sizes[1], MPI_DOUBLE, status.MPI_SOURCE, X_SOUTH, world);
    MPI_Send(y_south, sizes[1], MPI_DOUBLE, status.MPI_SOURCE, Y_SOUTH, world);
    delete i_south; delete x_south; delete y_south;
    //East stencil
    i_east = new int[sizes[2]];
    x_east = new double[sizes[2]];
    y_east = new double[sizes[2]];
    for(int i = 0; i < sizes[2]; i++){
        i_east[i] = s_index[East][i];
        x_east[i] = s_x[East][i];
        y_east[i] = s_y[East][i];
    }
    MPI_Send(i_east, sizes[2], MPI_INT, status.MPI_SOURCE, IND_EAST, world);
    MPI_Send(x_east, sizes[2], MPI_DOUBLE, status.MPI_SOURCE, X_EAST, world);
    MPI_Send(y_east, sizes[2], MPI_DOUBLE, status.MPI_SOURCE, Y_EAST, world);
    delete i_east; delete x_east; delete y_east;
    //West stencil
    i_west = new int[sizes[3]];
    x_west = new double[sizes[3]];
    y_west = new double[sizes[3]];
    for(int i = 0; i < sizes[3]; i++){
        i_west[i] = s_index[West][i];
        x_west[i] = s_x[West][i];
        y_west[i] = s_y[West][i];
    }
    MPI_Send(i_west, sizes[3], MPI_INT, status.MPI_SOURCE, IND_WEST, world);
    MPI_Send(x_west, sizes[3], MPI_DOUBLE, status.MPI_SOURCE, X_WEST, world);
    MPI_Send(y_west, sizes[3], MPI_DOUBLE, status.MPI_SOURCE, Y_WEST, world);
    delete i_west; delete x_west; delete y_west;
    
}
bool PDDSparseGM::Receive_Node(void){
    MPI_Comm_rank(workers, &workerid);
    MPI_Send(work_control, 2, MPI_INT, server, REQUEST_NODE, world);
    MPI_Recv(work_control, 2, MPI_INT, server, REPLY_NODE, world, MPI_STATUS_IGNORE);
    if(work_control[1] == 1){
        printf("Process %d is solving node %d .\n",workerid,work_control[0]);
        double pos[2];
        MPI_Recv(pos, 2, MPI_DOUBLE, server, NODE_POSITION, world, MPI_STATUS_IGNORE);
        position.resize(2);
        position[0] = pos[0]; position[1] = pos[1];
        printf("Located in [%.3f,%.3f]\n", position[0], position[1]);
        node_stencil = Recieve_Stencil_Data();
        node_stencil.Print(work_control[0]);
        return false;
    } else {
        printf("Process %d ended its work.\n",workerid);
        return true;
    }
}
Stencil PDDSparseGM::Recieve_Stencil_Data(void){
    double s_params[4];
    MPI_Recv(s_params, 4, MPI_DOUBLE, server, PARAMETERS, world, MPI_STATUS_IGNORE);
    int sizes[4];
    MPI_Recv(sizes, 4, MPI_INT, server, SIZES, world, MPI_STATUS_IGNORE);
    int *i_north, *i_south, *i_west, *i_east;
    double *x_north, *x_south, *x_west, *x_east, 
           *y_north, *y_south, *y_west, *y_east;
    std::map<direction, std::vector<int> > s_index;
    std::map<direction, std::vector<double> > s_x, s_y;
    i_north = new int[sizes[0]]; x_north = new double[sizes[0]]; y_north = new double[sizes[0]];
    MPI_Recv(i_north, sizes[0], MPI_INT, server, IND_NORTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_north, sizes[0], MPI_DOUBLE, server, X_NORTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_north, sizes[0], MPI_DOUBLE, server, Y_NORTH, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[0]; i++){
        s_index[North].push_back(i_north[i]);
        s_x[North].push_back(x_north[i]);
        s_y[North].push_back(y_north[i]);
    }
    delete i_north; delete x_north; delete y_north;
    i_south = new int[sizes[1]]; x_south = new double[sizes[1]]; y_south = new double[sizes[1]]; 
    MPI_Recv(i_south, sizes[1], MPI_INT, server, IND_SOUTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_south, sizes[1], MPI_DOUBLE, server, X_SOUTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_south, sizes[1], MPI_DOUBLE, server, Y_SOUTH, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[1]; i++){
        s_index[South].push_back(i_south[i]);
        s_x[South].push_back(x_south[i]);
        s_y[South].push_back(y_south[i]);
    }
    delete i_south; delete x_south; delete y_south;
    i_east = new int[sizes[2]]; x_east = new double[sizes[2]]; y_east = new double[sizes[2]];
    MPI_Recv(i_east, sizes[2], MPI_INT, server, IND_EAST, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_east, sizes[2], MPI_DOUBLE, server, X_EAST, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_east, sizes[2], MPI_DOUBLE, server, Y_EAST, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[2]; i++){
        s_index[East].push_back(i_east[i]);
        s_x[East].push_back(x_east[i]);
        s_y[East].push_back(y_east[i]);
    }
    delete i_east; delete x_east; delete y_east;
    i_west = new int[sizes[3]]; x_west = new double[sizes[3]]; y_west = new double[sizes[3]];
    MPI_Recv(i_west, sizes[3], MPI_INT, server, IND_WEST, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_west, sizes[3], MPI_DOUBLE, server, X_WEST, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_west, sizes[3], MPI_DOUBLE, server, Y_WEST, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[3]; i++){
        s_index[West].push_back(i_west[i]);
        s_x[West].push_back(x_west[i]);
        s_y[West].push_back(y_west[i]);
    }
    delete i_west; delete x_west; delete y_west;
    
    Stencil output;
    output.Init(s_index,s_x,s_y,s_params);
    return output;
}
void PDDSparseGM::Send_G_B(void){
    work_control[0] = (int) G.size();
    work_control[1] = (int) B.size();
    MPI_Send(work_control, 2, MPI_INT, server, REQUEST_NODE, world);
    //G_i is sent
    int *G_i_aux;
    G_i_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_i_aux[i] = G_i[i];
    G_i.resize(0);
    MPI_Send(G_i_aux, work_control[0], MPI_INT, server, TAG_Gi, world);
    delete G_i_aux;
    //G_j is sent
    int *G_j_aux;
    G_j_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_j_aux[i] = G_j[i];
    G_j.resize(0);
    MPI_Send(G_j_aux, work_control[0], MPI_INT, server, TAG_Gj, world);
    delete G_j_aux;
    //G_val is sent
    double *G_aux;
    G_aux = new double[work_control[0]];
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G[i];
    }
    G.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gval, world);
    delete G_aux;
    //B_i is sent
    int *B_i_aux;
    B_i_aux = new int[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_i_aux[i] = B_i[i];
    B_i.resize(0);
    MPI_Send(B_i_aux, work_control[1], MPI_INT, server, TAG_Bi, world);
    delete B_i_aux;
    //B is sent
    double *B_aux;
    B_aux = new double[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B[i];
    B.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bval, world);
    delete B_aux;
    std::cout << "Processor " << workerid << " ended its work \n";
    MPI_Comm_free(&workers);
}
void PDDSparseGM::Receive_G_B(void){
    //It ends worker process solving
    MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
    REQUEST_NODE, world, &status);
    work_control[0] = 0;
    work_control[1] = 0;
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    //The server accumulates the values of G and B
    MPI_Recv(work_control, 2, MPI_INT, status.MPI_SOURCE, REQUEST_NODE, world, &status);
    G.resize(work_control[0]);
    G_i.resize(work_control[0]);
    G_j.resize(work_control[0]);
    //G_i is received
    int *G_i_aux, *G_j_aux;
    double *G_aux;
    G_i_aux = new int[work_control[0]];
    MPI_Recv(G_i_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gi, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_i[i] = G_i_aux[i];
    delete G_i_aux;
    //G_j is received
    G_j_aux = new int[work_control[0]];
    MPI_Recv(G_j_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gj, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_j[i] = G_j_aux[i];
    delete G_j_aux;
    //G is received
    G_aux = new double[work_control[0]];
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gval, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    delete G_aux;
    //Triplets vector is fullfilled
    for(int j = 0; j < work_control[0]; j++){
        T_vec_G.push_back(T(G_i[j], G_j[j], G[j]));
    }
    G.resize(0);
    G_i.resize(0);
    G_j.resize(0);
    B.resize(work_control[1]);
    B_i.resize(work_control[1]);
    double *B_aux;
    int *B_i_aux;
    //B_i is received
    B_i_aux = new int[work_control[1]];
    MPI_Recv(B_i_aux, work_control[1], MPI_INT, status.MPI_SOURCE, TAG_Bi, world, &status);
    for(int i = 0; i< work_control[1]; i++) B_i[i] = B_i_aux[i];
    delete B_i_aux;
    //B is received
    B_aux = new double[work_control[1]];
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bval, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    delete B_aux;
    for(int j = 0; j < work_control[1]; j++){
        T_vec_B.push_back(T(B_i[j], 0, B[j]));
    }
}
Eigen::VectorXd PDDSparseGM::Node_Position(int node_index){
    Eigen::VectorXd output;
    output.resize(2);
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        if((*it).In(node_index)){
            output = (*it).Node_Position(node_index);
        }
    }
    return output;
}
void PDDSparseGM::Compute_Solution(BVP bvp){
    Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
    unsigned int size = T_vec_B.size();
    I_sparse.resize(size, size);
    I_sparse.setIdentity();
    G_sparse.resize(size, size);
    B_sparse.resize(size, 1);
    G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
    T_vec_G.resize(0);
    G_sparse+=I_sparse;
    I_sparse.resize(0,0);
    FILE *ofile;
    ofile = fopen("Output/Debug/G.csv","w");
    fprintf(ofile, "G_i,G_j,G_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
    T_vec_B.resize(0);
    ofile = fopen("Output/Debug/B.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    Eigen::VectorXd ud,Bd;
    Bd = Eigen::VectorXd(B_sparse);
    B_sparse.resize(0,0);
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_LU;
    solver_LU.compute(G_sparse);
    solver_LU.factorize(G_sparse);
    ud = solver_LU.solve(Bd);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_IT;
    solver_IT.compute(G_sparse);
    ud = solver_IT.solveWithGuess(Bd,ud);
    ofile = fopen("Output/solution.csv","w");
    fprintf(ofile,"Knot_index,x,y,sol_analytic,sol_PDDS,err,rerr\n");
    double sol, err, rerr;
    for(int i = 0; i< (int)ud.size(); i++){
        sol = bvp.u.Value(Node_Position(i),T_start);
        err = ud(i) - sol;
        rerr = fabs(err)/sol;
        fprintf(ofile,"%d,%f,%f,%f,%f,%f,%f\n",i,Node_Position(i)[0],Node_Position(i)[1],
        sol,ud(i),err,rerr);
    }
    fclose(ofile);
}
void PDDSparseGM::Print_Interface(void){
    if(myid==server){
    int cent = 0;
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        (*it).Print(cent);
        cent ++;
    }
    }
}
void PDDSparseGM::Print_Problem(void){
    printf("**********PDDSparse Solver*************\n");
    printf("***************Ver 1.0******************\n");
    printf("***********Jorge Moron Vidal************\n");
    printf("**********jmoron@math.uc3m.es***********\n");
    printf("h0 = %e \t Epsilon = %e \t T_start = %0.2f\n",h0, eps, T_start);
    printf("Number of trayectories per node %d\n",N);
    printf("Number of subdomains/interfaces: %d x %d \n", iN[0], iN[1]);
    printf("Number of nodes per interface: %d\n",nN[0]);
    printf("Position of South-West corner of the domain: (%.2f,%.2f)\n",SW(0), SW(1));
    printf("Position of North-East corner of the domain: (%.2f,%.2f)\n",NE(0), NE(1));
}
