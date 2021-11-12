//
//  ThirdPhase.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "ThirdPhase.hpp"
void Problem::THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage){
    //Create model
    mTHPMIP = IloModel (env);
    cplexmTHPMIP = IloCplex(mTHPMIP);
    cplexmTHPMIP.setParam(IloCplex::Param::Threads, 1);
    
    //Create decision variables and set decision variables values according to the scenario
    tcyvar = Create_tcyvar("r", Cycles2ndStage);
    tchvar = Create_tchvar("s", Chains2ndStage);
    
    //Create constraints
    IloRangeArray Disjoint(env);
    Disjoint = DisjointTHP(tcyvar, tchvar);
    
    //Solve formulation
    cplexmTHPMIP.solve();
    
    //Get THPMIP_Obj
        //if (THPMIP_Obj >= FPMIP_Obj){
        //  if (SPMIP_Obj < THPMIP_Obj){
        //      While (SPMIP_Obj < THPMIP_Obj){
        //          Add cols and constraints to SPMIP_Obj
        //      }
        //      if (SPMIP_Obj >= FPMIP_Obj){
        //          Robust Solution found
        //      }
        //      else{
        //          send scenario to FPMIP
        //      }
        //  }
        //}
        //else{
        //  send scenario to FPMIP
        //}
    
}
IloNumVarArray Problem::Create_tcyvar(const char* prefix, vector<Cycles>&Cycles2ndStage){
    //Get upper bound of cycles according to scenario
    map<int,bool>ub_tcyvar; // <cycle number, 0 or 1>
    ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
    
    IloNumVarArray var(env, Cycles2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        map<int,bool>::iterator find = ub_tcyvar.find(i);
        if (find != ub_tcyvar.end()){
            var[i] = IloNumVar(env, 0, 0, ILOINT);
            SetName(var[i], prefix, i + 1);
        }
        else{
            var[i] = IloNumVar(env, 0, 1, ILOINT);
            SetName(var[i], prefix, i + 1);
        }
    }
    return var;
}
map<int,bool> Problem::GetUB_tcyvar(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices){
    map<int,bool> map;
    for (auto it = FailedArcs.begin(); it != FailedArcs.end(); it++){
        for (auto it2 = ArcsinCyclesTHP[it->first].begin(); it2 != ArcsinCyclesTHP[it->first].end(); it2++){
            map[*(it2)] = true;
        }
    }
    for (auto it = FailedVertices.begin(); it != FailedVertices.end(); it++){
        for (auto it2 = CycleNodeTPH[it->first].begin(); it2 != CycleNodeTPH[it->first].end(); it2++){
            map[*(it2)] = true;
        }
    }
    return map;
}
IloNumVarArray Problem::Create_tchvar(const char* prefix, vector<Chain>&Chains2ndStage){
    //Get upper bound of cycles according to scenario
    map<int,bool>ub_tchvar; // <cycle number, 0 or 1>
    ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
    
    IloNumVarArray var(env, Chains2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        map<int,bool>::iterator find = ub_tchvar.find(i);
        if (find != ub_tchvar.end()){
            var[i] = IloNumVar(env, 0, 0, ILOINT);
            SetName(var[i], prefix, i + 1);
        }
        else{
            var[i] = IloNumVar(env, 0, 1, ILOINT);
            SetName(var[i], prefix, i + 1);
        }
    }
    return var;
}
map<int,bool> Problem::GetUB_tchvar(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices){
    map<int,bool> map;
    for (auto it = FailedArcs.begin(); it != FailedArcs.end(); it++){
        for (auto it2 = ArcsinChainsTHP[it->first].begin(); it2 != ArcsinChainsTHP[it->first].end(); it2++){
            map[*(it2)] = true;
        }
    }
    for (auto it = FailedVertices.begin(); it != FailedVertices.end(); it++){
        for (auto it2 = ChainNodeTPH[it->first].begin(); it2 != ChainNodeTPH[it->first].end(); it2++){
            map[*(it2)] = true;
        }
    }
    return map;
}
IloRangeArray Problem::DisjointTHP(IloNumVarArray& tcyvar, IloNumVarArray& tchvar){
    DisjointTHPArray = IloRangeArray(env);
    for (int i = 0; i < Nodes; i++){
        IloExpr expr (env,0);
        for (auto it = CycleNodeTPH[i].begin(); it!= CycleNodeTPH[i].end(); it++){
            expr += tcyvar[*(it)];
        }
        for (auto it = ChainNodeTPH[i].begin(); it!= ChainNodeTPH[i].end(); it++){
            expr += tchvar[*(it)];
        }
        DisjointTHPArray.add(IloRange(env, 1, expr, IloInfinity));
    }
    return DisjointTHPArray;
}
