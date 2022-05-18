//
//  main.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//


#include "main.hpp"


int main(int argc, const char * argv[]) { // cuatro parámetros ojo!!!
    if (argc != 12){
        cout <<"Not enough parameters. They must be 12." << endl;
        return -1;
    }
    
    //List of arguments
    string FolderName = argv[1];
    string FileName = argv[2];
    stringstream str; str << argv[3];
    IloInt CycleLength; str >> CycleLength;
    str.clear(); str << argv[4];
    IloInt ChainLength; str >> ChainLength;
    str.clear(); str << argv[5];
    IloInt VertexBudget; str >> VertexBudget;
    str.clear(); str << argv[6];
    IloInt ArcBudget; str >> ArcBudget;
    string RecoursePolicy = argv[7];
    string THP_Method = argv[8];
    string THP_Bound = argv[9];
    string WhereItisRun = argv[10];
    str.clear(); str << argv[11];
    IloNum TimeLimit; str >> TimeLimit;
    Problem P(FolderName, FileName, CycleLength, ChainLength, RecoursePolicy, THP_Method, THP_Bound, VertexBudget, ArcBudget, WhereItisRun, TimeLimit);
    P.ProgramStart = clock();
    
    cout << "Start reading" << endl;
    //clock_t tStart = clock();
    if (P.Reading() == 0) {
        cout << "Failed while reading input..." << endl;
        return -1;
    }
    
    P.ROBUST_KEP();
    cout << endl << "End" << endl;
    return 0;
}
