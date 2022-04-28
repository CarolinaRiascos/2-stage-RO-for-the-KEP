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
        direDonde = "../../../../../../0TestInstances/ROCarvalho2021/" + FolderName + "/";
    }
    else{
        direDonde = "../../0TestInstances/ROCarvalho2021/" + FolderName + "/";
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
void Problem::Print2ndStage(string status){
    //Print summary of results
    //"/Users/caroriascos/Documents/PhdProjectUofT/XcodeTests/Output/LagrangianBB/LagResults.txt"
    string OutputDire;
    if (WhereItisRun == "PC"){
        OutputDire = "../../../../../../Output/ROResults/TwoStageRO.txt";
    }else{
        OutputDire = "../Output/TwoStageROv4.txt";
    }
    file.open(OutputDire, fstream::in);
    if (this->file) {
        this->file.close();
        this->file.open(OutputDire, fstream::app);
    }else{
        file.open(OutputDire, fstream::out);
        file << "Instance" << '\t' << "PDP" << '\t' << "NDD" << '\t' << "K" << '\t' << "L" << '\t' << "VertexBudget" << '\t' << "ArcBudget" << '\t' << "Method" << '\t' << "Policy" << '\t' << "status" << '\t' << "Ite1stS" << '\t' << "TotalTime1stS" << '\t' << "TotalIte2ndS" << '\t' << "TotalTime2ndS" << '\t' << "TotalTimeMP2ndS" << '\t' << "TotalTimeHeu" << '\t' << "TotalIteHeuTrue" <<  '\t' << "TotalTimeRecoCG" << '\t' << "TotalIteCGTrue" << '\t' << "TotalTimeRecoMIP" << '\t' << "RO_Objective" << '\t' << "TotalIteOptP" << '\t' << "TotalIteOptPIte1stis1" << '\t' << "TotalTimeOptP" << '\t' << "nOutInfeas" << '\t' << "nOutBound" << endl;
    }
    
    file << FileName << '\t' << Pairs << '\t' << NDDs << '\t' << CycleLength << '\t' << ChainLength << '\t' << MaxVertexFailures <<'\t' << MaxArcFailures << '\t' << THP_Method << '\t' << RecoursePolicy << '\t' << status << '\t' << Ite1stStage << '\t' << tTotal1stS << '\t' << GlobalIte2ndStage << '\t' << tTotal2ndS << '\t' << tTotalMP2ndPH << '\t' << tTotalHeu << '\t' << runHeuristicstrue <<  '\t' << tTotalRecoCG << '\t' << runCGtrue << '\t' << tTotalRecoMIP << '\t' << FPMIP_Obj << '\t' << IteOptP << '\t' << IteOptPIte1stis1 << '\t' << tTotalOptP << '\t' << OutforInfeas << '\t' << OutforBound << endl;
    
    file.close();
    exit(0);
}
