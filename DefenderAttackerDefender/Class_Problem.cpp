//
//  Class_Problem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Class_Problem.hpp"
Problem::Problem(string _FolderName, string _FileName, IloInt _cycleLength, IloInt _chainLength, string _RecoursePolicy){
    FolderName = _FolderName;
    FileName = _FileName;
    CycleLength = _cycleLength;
    ChainLength = _chainLength;
    RecoursePolicy = _RecoursePolicy;
}
