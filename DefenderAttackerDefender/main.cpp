//
//  main.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//


#include "main.hpp"


int main(int argc, const char * argv[]) { // cuatro parámetros ojo!!!
    if (argc != 6){
        cout <<"Not enough parameters. They must be 6." << endl;
        return -1;
    }
    
    
    //List of arguments
    string FolderName = argv[1];
    string FileName = argv[2];
    stringstream str; str << argv[3];
    IloInt CycleLength; str >> CycleLength;
    str.clear(); str << argv[4];
    IloInt ChainLength; str >> ChainLength;
    string RecoursePolicy = argv[5];
    Problem P(FolderName, FileName, CycleLength, ChainLength, RecoursePolicy);
    
    cout << "Start reading" << endl;
    //clock_t tStart = clock();
    if (P.Reading() == 0) {
        cout << "Failed while reading input..." << endl;
        return -1;
    }
    
    
    P.M_PICEF();
    cout << endl << "End" << endl;
    return 0;
}
