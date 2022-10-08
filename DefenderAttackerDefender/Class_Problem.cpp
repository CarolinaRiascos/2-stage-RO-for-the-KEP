//
//  Class_Problem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Class_Problem.hpp"
Problem::Problem(string _FilePath,  IloInt _cycleLength, IloInt _chainLength, string _RecoursePolicy, string _THP_Method, string _THP_Bound, IloInt _VertexBudget, IloInt _ArcBudget, string _OutputPath, IloNum _TimeLimit){
    FilePath = _FilePath;
    CycleLength = _cycleLength;
    ChainLength = _chainLength;
    RecoursePolicy = _RecoursePolicy;
    THP_Method = _THP_Method;
    THP_Bound = _THP_Bound;
    MaxVertexFailures = _VertexBudget;
    MaxArcFailures = _ArcBudget;
    OutputPath = _OutputPath;
    TimeLimit = _TimeLimit;
}
void Problem::SetName(IloNumVar& var, const char* prefix, IloInt i){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i);
    const char* varName = name.c_str();
    var.setName(varName);
}
void Problem::SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
