//
//  Reading.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Reading.hpp"

int Problem::Reading() {

    ifstream inFile(FilePath, ifstream::in);
    if (!inFile) {
        cout << endl << "Instance's files not found. " << FilePath << endl;
        return 0;
    }
    //Get FileName
    string new_name = FilePath; long pos_str = -1; long size_str = FilePath.size();
    while(pos_str <= size_str){
        new_name = new_name.substr(pos_str + 1);
        pos_str = new_name.find("/");
        if (pos_str < 0) break;
    }
    FileName = new_name;
    
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
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                PredMap[AdjacencyList[i][j] - 1].push_back(make_pair(i, j));
        }
    }
    
    return 1;
}
void Problem::Print2ndStage(string status, vector<IndexGrandSubSol>&SolFirstStage){
    //Print summary of results
    //"/Users/caroriascos/Documents/PhdProjectUofT/XcodeTests/Output/LagrangianBB/LagResults.txt"
    cout << to_string(LeftTime) << ": " <<  status << endl;
    
    file.open(OutputPath, fstream::in);
    if (this->file) {
        this->file.close();
        this->file.open(OutputPath, fstream::app);
    }else{
        file.open(OutputPath, fstream::out);
        file << "Instance" << '\t' << "PDP" << '\t' << "NDD" << '\t' << "K" << '\t' << "L" << '\t' << "VertexBudget" << '\t' << "ArcBudget" << '\t' << "Method" << '\t' << "LB" << '\t' << "Policy" << '\t' << "status" << '\t' << "Ite1stS" << '\t' << "TotalTime1stS" << '\t' << "TotalIte2ndS" << '\t' << "TotalTime2ndS" << '\t' << "TotalTimeMP2ndS" << '\t' << "TotalTimeHeu" << '\t' << "TotalIteHeuTrue" <<  '\t' << "TotalTimeRecoCG" << '\t' << "TotalIteCGTrue" << '\t' << "TotalTimeRecoMIP" << '\t' << "RO_Objective" << '\t' << "TotalIteOptP" << '\t' << "TotalIteOptPIte1stis1" << '\t' << "TotalTimeOptP" << '\t' << "nOutInfeas" << '\t' << "nOutBound" << '\t' << "TimeFindingCyCh" << '\t' << "NumDominatedS" << endl;
    }
    
    file << FileName << '\t' << Pairs << '\t' << NDDs << '\t' << CycleLength << '\t' << ChainLength << '\t' << MaxVertexFailures <<'\t' << MaxArcFailures << '\t' << THP_Method << '\t' << THP_Bound << '\t' << RecoursePolicy << '\t' << status << '\t' << Ite1stStage << '\t' << tTotal1stS << '\t' << GlobalIte2ndStage << '\t' << tTotal2ndS << '\t' << tTotalMP2ndPH << '\t' << tTotalHeu << '\t' << runHeuristicstrue <<  '\t' << tTotalRecoCG << '\t' << runCGtrue << '\t' << tTotalRecoMIP << '\t' << FPMIP_Obj << '\t' << IteOptP << '\t' << IteOptPIte1stis1 << '\t' << tTotalOptP << '\t' << OutforInfeas << '\t' << OutforBound << '\t' << tTotalFindingCyCh << '\t' << NumDominatedS << endl;
    
    cout << "First-stage solution:";
    int ncy = 0; int nch = 0;
    for (int i = 0; i < SolFirstStage.size(); i++){
        cout << endl;
        if (SolFirstStage[i].get_cc()[0] < Pairs){
            ncy++;
            cout << "Cycle " << to_string(ncy) << ":" << "\t";
        }
        else{
            nch++;
            cout << "Chain " << to_string(nch) << ":" << "\t";
        }
        for (int j = 0; j < SolFirstStage[i].get_cc().size(); j++){
            cout << SolFirstStage[i].get_cc()[j] + 1 << "\t";
        }
    }
//    cout << endl << "Recourse solution: " << endl;
//    ncy = 0; nch = 0;
//    for (int i = 0; i < RecoMatching.size(); i++){
//        cout << endl;
//        if (RecoMatching[i][0] < Pairs){
//            ncy++;
//            cout << "Cycle " << to_string(ncy) << ":" << "\t";
//        }
//        else{
//            nch++;
//            cout << "Chain " << to_string(nch) << ":" << "\t";
//        }
//        for (int j = 0; j < RecoMatching[i].size(); j++){
//            cout << RecoMatching[i][j] + 1 << "\t";
//        }
//    }
    
    cout << endl << to_string(LeftTime) + ": Results printed" << endl;
    
    file.close();
    exit(0);
}
