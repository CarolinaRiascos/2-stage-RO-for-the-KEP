//
//  GrandSubProbMasterPICEF.hpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2022-05-05.
//  Copyright © 2022 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef GrandSubProbMasterPICEF_hpp
#define GrandSubProbMasterPICEF_hpp

#include <stdio.h>
#include "Class_Problem.hpp"
NumVar4D Create4DBinTHP(IloEnv& env, vector<IndexGrandSubSol>& ChainsTHP, int L, string name);
bool wasSelected(vector<int>&v, int a);
IloExpr GetObjTPH_PICEF(IloEnv& env, vector<Cycles>&Cycles2ndStage, vector<int>&SelectedVertices, IloNumVarArray& var1D, NumVar2D& var2D, NumVar3D& var3D, IloNumArray2& AdjaList, int Pairs, int L, map<int,bool>& ub_tcyvar);
#endif /* GrandSubProbMasterPICEF_hpp */
