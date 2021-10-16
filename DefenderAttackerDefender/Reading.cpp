//
//  Reading.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Reading.hpp"

int Problem::Reading() {

    string direDonde = "/Users/caroriascos/Documents/PhdProjectUofT/XcodeTests/0TestInstances/InstancesFormat/" + FolderName + "/";
    //
    //"/home/criascos/codes/LagBB/Instances/" + FolderName + "/";

    ifstream inFile(direDonde + FileName, ifstream::in);
    if (!inFile) {
        cout << endl << "Intance's files not found. " << FileName << endl;
        return 0;
    }
    AdjacencyList = IloNumArray2(env);
    WeightMatrix = IloNumArray2(env);
    //Leer Argumentos
    inFile >> Nodes >> NDDs >> Pairs >> NumArcs >> AdjacencyList >> WeightMatrix;


    //Build weights matrix
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            Weights[make_pair(i,AdjacencyList[i][j] - 1)] = WeightMatrix[i][j];
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
