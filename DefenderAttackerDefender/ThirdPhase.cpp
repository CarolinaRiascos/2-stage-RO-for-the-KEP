//
//  ThirdPhase.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "ThirdPhase.hpp"
void Problem::THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<int>&ListSelVertices){
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
    
    //Retrieve Solution and check whether the cycles are contained in Cycles2ndTo3rd and Chains2ndTo3rd
    IloNumArray tcysol (env, tcyvar.getSize());
    cplexmTHPMIP.getValues(tcysol,tcyvar);
    IloNumArray tchsol (env, tcyvar.getSize());
    cplexmTHPMIP.getValues(tchsol,tchvar);
    AddNewCols3rdTo2nd (tcysol, tchsol, Cycles2ndTo3rd, Chains2ndTo3rd, Cycles3rdTo2nd, Chains3rdTo2nd, Cyclenewrow2ndPH, Chainnewrow2ndPH, tcysol3rd, tchsol3rd, Cycles2ndStage, Chains2ndStage);
    
    for (int i = 0; i < tcysol3rd.size(); i++){
        //OldActiveCCSubP_CY
        vector<int> OldActiveCY;
        OldActiveCY = ModifyOldActiveCCSubP_CY(tcysol3rd[i], Cycles2ndTo3rd, Cyclenewrow2ndPH, Cycles2ndStage);
        //OldSelVert2ndPH
        vector<int> OldSelVert2ndPH;
        OldSelVert2ndPH = ModifyOldSelectVex_CY(tcysol3rd[i], ListSelVertices, Cycles2ndStage);

        //Create new cycle column
        IloNumColumn col(env);
        //Create new variable
        cyvar.add(IloNumVar(col));
        string name = "x." + to_string(Cyclenewrow2ndPH + i);
        const char* varName = name.c_str();
        cyvar[cyvar.getSize() - 1].setName(varName);
        cyvar[cyvar.getSize() - 1].setBounds(0, 1);
        
        //Update Beta constraint
        col+= ConsBeta[0](-1);
        //UpdateSelVert2ndPH
        for (int i = 0; i < OldSelVert2ndPH.size(); i++){
            SelVert2ndPH[OldSelVert2ndPH[i]](-1);
        }
        //UpdateActiveCCSubP_CY
        for (int i = 0; i < OldActiveCY.size(); i++){
            ActiveCCSubP_CY[OldActiveCY[i]](-1);
        }
        //NewActiveCCSubP_CY
        GetNewIloRangeCY3rd(tcysol3rd[i], Cycles2ndStage);
        GrandSubProb.add(ActiveCCSubP_CY[ActiveCCSubP_CY.getSize() - 1]);
        
    }
    for (int i = 0; i < tchsol3rd.size(); i++){
        //OldActiveCCSubP_CH
        vector<int> OldActiveCH;
        OldActiveCH = ModifyOldActiveCCSubP_CH(tchsol3rd[i], Cycles2ndTo3rd, Chainnewrow2ndPH, Chains2ndStage);
        //OldSelVert2ndPH
        vector<int> OldSelVert2ndPH;
        OldSelVert2ndPH = ModifyOldSelectVex_CH(tchsol3rd[i], ListSelVertices, Chains2ndStage);
        //Create new chain column
        IloNumColumn col(env);
        //Create new variable
        chvar.add(IloNumVar(col));
        string name = "y." + to_string(Chainnewrow2ndPH + i);
        const char* varName = name.c_str();
        chvar[chvar.getSize() - 1].setName(varName);
        chvar[chvar.getSize() - 1].setBounds(0, 1);
        
        //Update Beta constraint
        col+= ConsBeta[0](-1);
        //UpdateSelVert2ndPH
        for (int i = 0; i < OldSelVert2ndPH.size(); i++){
            SelVert2ndPH[OldSelVert2ndPH[i]](-1);
        }
        //UpdateActiveCCSubP_CH
        for (int i = 0; i < OldActiveCH.size(); i++){
            ActiveCCSubP_CH[OldActiveCH[i]](-1);
        }
        //NewActiveCCSubP_CH
        GetNewIloRangeCH3rd(tcysol3rd[i], Chains2ndStage);
        GrandSubProb.add(ActiveCCSubP_CH[ActiveCCSubP_CH.getSize() - 1]);
    }
    //Add NoGoodCut
    
    
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
void Problem::AddNewCols3rdTo2nd (IloNumArray tcysol, IloNumArray tchsol, map<int,int>& Cycles2ndTo3rd, map<int,int>& Chains2ndTo3rd, map<int,int>& Cycles3rdTo2nd, map<int,int>& Chains3rdTo2nd, int& Cyclenewrow2ndPH, int& Chainnewrow2ndPH, vector<int>& tcysol2nd, vector<int>& tchsol2nd, vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage){
    Cyclenewrow2ndPH = int(Cycles2ndTo3rd.size());
    Chainnewrow2ndPH = int(Chains2ndTo3rd.size());
    tcysol2nd.clear();
    tchsol2nd.clear();
    for (int i = 0; i < tcysol.getSize(); i++){
        if (tcysol[i] > 0.9){
            tcysol2nd.push_back(i);
            map<int,int>::iterator it = Cycles3rdTo2nd.find(i);
            if (it == Cycles3rdTo2nd.end()){
                Cycles2ndTo3rd[int(Cycles2ndTo3rd.size() - 1)] = i;
                Cycles3rdTo2nd[i] = int(Cycles2ndTo3rd.size() - 1);
                //Update CycleNodeSPH
                for (int j = 0; j < Cycles2ndStage[i].get_c().size(); j++){
                    CycleNodeSPH[Cycles2ndStage[i].get_c()[j]].push_back(int(Cycles2ndTo3rd.size() - 1));
                }
            }
        }
    }
    for (int i = 0; i < tchsol.getSize(); i++){
        if (tchsol[i] > 0.9){
            tchsol2nd.push_back(i);
            map<int,int>::iterator it = Chains3rdTo2nd.find(i);
            if (it == Chains3rdTo2nd.end()){
                Chains2ndTo3rd[int(Chains2ndTo3rd.size() - 1)] = i;
                Chains3rdTo2nd[i] = int(Chains2ndTo3rd.size() - 1);
                //Update ChainNodeSPH
                for (int j = 0; j < Chains2ndStage[i].Vnodes.size(); j++){
                    CycleNodeSPH[Chains2ndStage[i].Vnodes[j].vertex].push_back(int(Chains2ndTo3rd.size() - 1));
                }
            }
        }
    }
    
}
vector<int>Problem::ModifyOldActiveCCSubP_CY(int tOnecysol, map<int,int>& Cycles2ndTo3rd, int& Cyclenewrow2ndPH, vector<Cycles>&Cycles2ndStage){
    vector<int>UpdateActiveCCSubP_CY;
    for (int i = 0; i < Cyclenewrow2ndPH; i++){
        bool found = false;
        for (int j = 0; j < Cycles2ndStage[Cycles2ndTo3rd[i]].get_c().size(); j++){
            int u = Cycles2ndStage[Cycles2ndTo3rd[i]].get_c()[j];
            for (int h = 0; h < Cycles2ndStage[tOnecysol].get_c().size(); h++){
                int v = Cycles2ndStage[tOnecysol].get_c()[h];
                if (u == v){
                    UpdateActiveCCSubP_CY.push_back(i);
                    found = true;
                    break;
                }
            }
            if (found == true) break;
        }
    }
    return UpdateActiveCCSubP_CY;
}
vector<int>Problem::ModifyOldActiveCCSubP_CH(int tOnechsol, map<int,int>& Chains2ndTo3rd, int& Chainnewrow2ndPH, vector<Chain>&Chains2ndStage){
    vector<int>UpdateActiveCCSubP_CH;
    for (int i = 0; i < Chainnewrow2ndPH; i++){
        bool found = false;
        for (int j = 0; j < Chains2ndStage[Chains2ndTo3rd[i]].Vnodes.size(); j++){
            int u = Chains2ndStage[Chains2ndTo3rd[i]].Vnodes[j].vertex;
            for (int h = 0; h < Chains2ndStage[tOnechsol].Vnodes.size(); h++){
                int v = Chains2ndStage[tOnechsol].Vnodes[h].vertex;
                if (u == v){
                    UpdateActiveCCSubP_CH.push_back(i);
                    found = true;
                    break;
                }
            }
            if (found == true) break;
        }
    }
    
    return UpdateActiveCCSubP_CH;
}
vector<int>Problem::ModifyOldSelectVex_CY(int tOnecysol2nd, vector<int>ListSelVertices, vector<Cycles>&Cycles2ndStage){
    vector<int>UpdateOldSelectVex_CY;
    
    for (int i = 0; i < ListSelVertices.size(); i++){
        for (int j = 0; j < Cycles2ndStage[tOnecysol2nd].get_c().size(); j++){
            if (Cycles2ndStage[tOnecysol2nd].get_c()[j] == ListSelVertices[i]){
                UpdateOldSelectVex_CY.push_back(i);
            }
        }
    }
    
    return UpdateOldSelectVex_CY;
}
vector<int>Problem::ModifyOldSelectVex_CH(int tOnechsol2nd, vector<int>ListSelVertices, vector<Chain>&Chains2ndStage){
    vector<int>UpdateOldSelectVex_CH;
    
    for (int i = 0; i < ListSelVertices.size(); i++){
        for (int j = 0; j < Chains2ndStage[tOnechsol2nd].Vnodes.size(); j++){
            if (Chains2ndStage[tOnechsol2nd].Vnodes[j].vertex == ListSelVertices[i]){
                UpdateOldSelectVex_CH.push_back(i);
            }
        }
    }
    
    return UpdateOldSelectVex_CH;
}
IloRange Problem::GetNewIloRangeCY3rd(int tOnecysol3rd, vector<Cycles>&Cycles2ndStage){
    
    IloExpr ExprArcsVtx (env, 0);
    IloExpr AllccVars (env, 0);
    map<int,bool>forcycles;
    map<int,bool>forchains;
    int idx = tOnecysol3rd;
    for (int j = 0; j < Cycles2ndStage[idx].get_c().size(); j++){
        ExprArcsVtx+= vertex[Cycles2ndStage[idx].get_c()[j]];
        int u = Cycles2ndStage[idx].get_c()[j];
        if (j == Cycles2ndStage[idx].get_c().size() - 1){
            int v = Cycles2ndStage[idx].get_c()[0];
            int s = mapArcs[make_pair(u,v)];
            ExprArcsVtx+= arc[u][s];
        }
        else{
            int v = Cycles2ndStage[idx].get_c()[j + 1];
            int s = mapArcs[make_pair(u,v)];
            ExprArcsVtx+= arc[u][s];
        }
        for (int k = 0; k < CycleNodeSPH[u].size(); k++){
            if (CycleNodeSPH[u][k] != Cycles3rdTo2nd[idx]) forcycles[CycleNodeSPH[u][k]];
            //check whether CycleNodeSPH[Cycles2ndStage[idx].get_c()[j]][k] exixts in Cycles2ndTo3rd
        }
        //cout << Cycles2ndStage[idx].get_c()[j] << endl;
        for (int k = 0; k < ChainNodeSPH[u].size(); k++){
            forchains[ChainNodeSPH[u][k]];
        }
    }
    for (auto it = forcycles.begin(); it != forcycles.end(); it++){
        AllccVars+= cyvar[it->first];
    }
    for (auto it = forchains.begin(); it != forchains.end(); it++){
        AllccVars+= chvar[it->first];
    }
    string name = "Active_cyvar2nd." + to_string(Cycles3rdTo2nd[idx]);
    const char* cName = name.c_str();
    //cout << cyvar[i].getName()  << ">=" <<  1 - ExprArcsVtx - AllccVars << endl;
    ActiveCCSubP_CY.add(IloRange(env, 1, cyvar[Cycles3rdTo2nd[idx]] + ExprArcsVtx + AllccVars, IloInfinity, cName));
    
    return NewIloRangeCY3rd;
}
IloRange Problem::GetNewIloRangeCH3rd(int tOnecysol3rd, vector<Chain>&Chains2ndStage){
    IloExpr Expr2ArcsVtx (env, 0);
    IloExpr All2ccVars (env, 0);
    map<int,bool>forcycles;
    map<int,bool>forchains;
    int idx = tOnecysol3rd;
    for (int j = 0; j < Chains2ndStage[idx].Vnodes.size(); j++){
        int u = Chains2ndStage[idx].Vnodes[j].vertex;
        Expr2ArcsVtx+= vertex[u];
        if (j <= Chains2ndStage[idx].Vnodes.size() - 2){
            int v = Chains2ndStage[idx].Vnodes[j + 1].vertex;
            int s = mapArcs[make_pair(u,v)];
            Expr2ArcsVtx+= arc[u][s];
        }
        for (int k = 0; k < CycleNodeSPH[u].size(); k++){
            forcycles[CycleNodeSPH[u][k]];
        }
        for (int k = 0; k < ChainNodeSPH[u].size(); k++){
            if (ChainNodeSPH[u][k] != Cycles3rdTo2nd[idx]) forchains[ChainNodeSPH[u][k]];
        }
    }
    for (auto it = forcycles.begin(); it != forcycles.end(); it++){
        All2ccVars+= cyvar[it->first];
    }
    for (auto it = forchains.begin(); it != forchains.end(); it++){
        All2ccVars+= chvar[it->first];
    }
    string name = "Active_chvar2nd." + to_string(Cycles3rdTo2nd[idx]);
    const char* cName = name.c_str();
    //cout << cyvar[i].getName()  << ">=" <<  1 - Expr2ArcsVtx - All2ccVars << endl;
    ActiveCCSubP_CH.add(IloRange(env, 1, chvar[Cycles3rdTo2nd[idx]] + Expr2ArcsVtx + All2ccVars, IloInfinity, cName));
    
    return NewIloRangeCH3rd;
}
