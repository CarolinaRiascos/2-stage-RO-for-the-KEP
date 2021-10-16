//
//  M-PICEF.hpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-09-22.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef M_PICEF_hpp
#define M_PICEF_hpp

#include <stdio.h>
#include "Class_Problem.hpp"


ILOSTLBEGIN

void SetName(IloNumVar& var, const char* prefix, IloInt i);
void SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j);
void SetName3(IloNumVar& var, const char* prefix, IloInt e, IloInt i, IloInt j);
bool inFVS (vector<int> v, int vertex);
int posinFVS(vector<int>fvs, int lookfor);
bool checkIfArcFeasible(map<pair<int,int>,vector<int>>& ArcFvs, pair<int,int>pair);



#endif /* M_PICEF_hpp */
