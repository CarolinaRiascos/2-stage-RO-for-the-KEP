//
//  Class_Problem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Class_Problem.hpp"
Problem::Problem(string _FolderName, string _FileName, IloInt _cycleLength, IloInt _chainLength, string _RecoursePolicy, string _THP_Method, IloInt _VertexBudget, IloInt _ArcBudget, string _WhereItisRun, IloNum _TimeLimit){
    FolderName = _FolderName;
    FileName = _FileName;
    CycleLength = _cycleLength;
    ChainLength = _chainLength;
    RecoursePolicy = _RecoursePolicy;
    THP_Method = _THP_Method;
    MaxVertexFailures = _VertexBudget;
    MaxArcFailures = _ArcBudget;
    WhereItisRun = _WhereItisRun;
    TimeLimit = _TimeLimit;
}
