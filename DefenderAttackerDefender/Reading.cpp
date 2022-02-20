//
//  Reading.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Reading.hpp"

int Problem::Reading() {
    string direDonde;
    if (WhereItisRun == "PC"){
        direDonde = "../../../../../../0TestInstances/PrefLib/" + FolderName + "/";
    }
    else{
        direDonde = "../../0TestInstances/PrefLib/" + FolderName + "/";
    }
    
    //
    //"/home/criascos/codes/LagBB/Instances/" + FolderName + "/";

    ifstream inFile(direDonde + FileName, ifstream::in);
    if (!inFile) {
        cout << endl << "Intance's files not found. " << FileName << endl;
        return 0;
    }
    AdjacencyList = IloNumArray2(env);
    WeightMatrix = IloNumArray2(env);
    PRAList = IloNumArray(env);
    //Leer Argumentos
    inFile >> Nodes >> NDDs >> Pairs >> NumArcs >> AdjacencyList >> WeightMatrix >> PRAList;


    //Build weights matrix
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            if (WeightMatrix.getSize() == 1){
                Weights[make_pair(i,AdjacencyList[i][j] - 1)] = 1;
            }
            else{
                Weights[make_pair(i,AdjacencyList[i][j] - 1)] = WeightMatrix[i][j];
            }
        }
    }
    
    //Build Predeccessors List
    PredList = vector<vector<int>>(AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                //PredList[i].push_back(j);
                PredMap[AdjacencyList[i][j] - 1].push_back(make_pair(i, j));
        }
    }
    
    return 1;
}
void Problem::Print2ndStage(){
    //Print summary of results
    //"/Users/caroriascos/Documents/PhdProjectUofT/XcodeTests/Output/LagrangianBB/LagResults.txt"
    string OutputDire;
    if (WhereItisRun == "PC"){
        OutputDire = "../../../../../../Output/ROResults/2ndStageRO.txt";
    }else{
        OutputDire = "../Output/2ndStageRO.txt";
    }
    file.open(OutputDire, fstream::in);
    if (this->file) {
        this->file.close();
        this->file.open(OutputDire, fstream::app);
    }else{
        file.open(OutputDire, fstream::out);
        file << "Instance" << '\t' << "PDP" << '\t' << "NDD" << '\t' << "K" << '\t' << "L" << '\t' << "VertexBudget" << '\t' << "ArcBudget" << '\t' <<"Ite" << '\t' << "Time(s)" << '\t' << "Method" << '\t' << "Policy" << '\t' << "SPMIP_Obj" << endl;
    }
    
    file << FileName << '\t' << Pairs << '\t' << NDDs << '\t' << CycleLength << '\t' << ChainLength << '\t' << MaxVertexFailures <<'\t' << MaxArcFailures <<'\t' << Ite2ndS << '\t' << tEnd2ndS << '\t' << THP_Method << '\t' << RecoursePolicy << '\t' << SPMIP_Obj << endl;
    
    file.close();
}
