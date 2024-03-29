//
//  main.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//


#include "main.hpp"


int main(int argc, const char * argv[]) { // cuatro parámetros ojo!!!
    if (argc != 11){
        cout <<"Not enough parameters. They must be 11." << endl;
        return -1;
    }
    
    //List of arguments
    string FilePath = argv[1];
    //string FileName = argv[2];
    stringstream str; str << argv[2];
    IloInt CycleLength; str >> CycleLength;
    str.clear(); str << argv[3];
    IloInt ChainLength; str >> ChainLength;
    str.clear(); str << argv[4];
    IloInt VertexBudget; str >> VertexBudget;
    str.clear(); str << argv[5];
    IloInt ArcBudget; str >> ArcBudget;
    string RecoursePolicy = argv[6];
    string THP_Method = argv[7];
    string THP_Bound = argv[8];
    string OutputPath = argv[9];
    str.clear(); str << argv[10];
    IloNum TimeLimit; str >> TimeLimit;
    Problem P(FilePath, CycleLength, ChainLength, RecoursePolicy, THP_Method, THP_Bound, VertexBudget, ArcBudget, OutputPath, TimeLimit);
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
