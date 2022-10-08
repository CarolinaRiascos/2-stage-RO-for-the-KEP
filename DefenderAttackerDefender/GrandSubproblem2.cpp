//
//  GrandSubproblem2.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-12-03.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandSubproblem2.hpp"
void Problem::GrandSubProbMaster2(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<IndexGrandSubSol>&SolFirstStage){
    // Create model
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    GrandSubProb = IloModel(env);
    cplexGrandSubP = IloCplex(GrandSubProb);
    cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexGrandSubP.setParam(IloCplex::Param::Threads, 1);
    cplexGrandSubP.setOut(env.getNullStream());
    tStart2ndS = clock();
    
    KEPSols2ndStage.clear();
    SampleCols2ndStage2(Chains2ndStage, Cycles2ndStage, SolFirstStage);
    
    //Get selected vertices
    vector<int>ListSelVertices = GetSelPairs(SolFirstStage);
    
    //Create cycle variables
    cyvar = IloNumVarArray(env, Cycles2ndTo3rd.size(), 0, 1, ILOINT);
    for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
       SetName1Index(cyvar[i], "x", i + 1);
       //cout << r[i].getName() << endl;
    }
    
    //Create chain variables
    chvar = IloNumVarArray(env, Chains2ndTo3rd.size(), 0, 1, ILOINT);
    for (int i = 0; i < Chains2ndTo3rd.size(); i++){
       SetName1Index(chvar[i], "y", i + 1);
       //cout << r[i].getName() << endl;
    }
    
    //Create arc variables
    mapArcs.clear();
    arc = NumVar2D(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        arc[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, 1, ILOINT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            mapArcs[make_pair(i, AdjacencyList[i][j] - 1)] = j;
            SetName2(arc[i][j], "arc", i + 1, AdjacencyList[i][j]);
        }
    }
    
    //Create vertex variables
    vertex = IloNumVarArray(env, Nodes, 0, 1, ILOINT);
    for (int i = 0; i < Nodes; i++){
        SetName1Index( vertex[i], "vertex",i + 1);
    }
    
    //Create bounding variable
    Beta = IloNumVar(env);
    SetName1Index(Beta, "Beta", 0);
    
    //Create Bounding Constraint
    vBoundConstraint = IloRangeArray(env);
    Const11b(KEPSols2ndStage);
    GrandSubProb.add(vBoundConstraint[vBoundConstraint.getSize() - 1]);
    exprBound.end();
    
    //Create Active Cols Constraint
    vFailedMatches = IloRangeArray(env);
    Const11c(KEPSols2ndStage);
//    GrandSubProb.add(vFailedMatches);
//    exprVxtArcsCH.end();
//    exprVxtArcsCY.end();
    
    //Create arcs and vertex constraints
    IloExpr sumVertices (env, 0);
    for (int i = 0; i < Nodes; i++){
        sumVertices+=  vertex[i];
    }
    string name;
    const char* cName;
    name = "VtxSum";
    cName = name.c_str();
    GrandSubProb.add(IloRange(env, -IloInfinity, sumVertices, MaxVertexFailures, cName));
    
    IloExpr sumArcs (env, 0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumArcs+= arc[i][j];
        }
    }
    name = "ArcSum";
    cName = name.c_str();
    GrandSubProb.add(IloRange(env, -IloInfinity, sumArcs, MaxArcFailures, cName));
    
    //Add objective
    name = "Obj_GrandSubP";
    ObjGrandSubP = IloObjective(env, Beta, IloObjective::Minimize, name.c_str());
    GrandSubProb.add(ObjGrandSubP);
    
    //cplexGrandSubP.exportModel("GrandSubP2.lp");
    cplexGrandSubP.solve();
    if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
        cout << "S.O.S. This should not happen." << endl;
    }
    else{
        SPMIP_Obj = cplexGrandSubP.getValue(Beta);
        //env.out() << "Objective: " << SPMIP_Obj << endl;//cplexGrandSubP.getValue(Excess)
        
        //Retrieve solution
        IloNumArray cyvar_sol(env, Cycles2ndTo3rd.size());
        IloNumArray chvar_sol(env, Chains2ndTo3rd.size());
        cplexGrandSubP.getValues(cyvar_sol,cyvar);
        cplexGrandSubP.getValues(chvar_sol,chvar);
        
//        for (int i = 0; i < cyvar_sol.getSize(); i++){
//            if (cyvar_sol[i] > 0.9) cout << cyvar[i].getName() << endl;
//        }
//        for (int i = 0; i < chvar_sol.getSize(); i++){
//            if (chvar_sol[i] > 0.9) cout << chvar[i].getName() << endl;
//        }
        //cout << "Beta: " << cplexGrandSubP.getValue(Beta) << endl;
        
        vertex_sol = IloNumArray(env, Nodes);
        cplexGrandSubP.getValues(vertex_sol,vertex);
        
        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            cplexGrandSubP.getValues(arc_sol[f],arc[f]);
        }
        THPMIP(Cycles2ndStage, Chains2ndStage, SolFirstStage);
    }

}
void Problem::Const11b(vector<KEPSol>&KEPSols2ndStage){
    
    exprBound = IloExpr(env,0);
    for (int i = 0; i < KEPSols2ndStage.back().cycles.size(); i++){
        int w = Cycles2ndStage[Cycles2ndTo3rd[KEPSols2ndStage.back().cycles[i]]].get_Many();
        exprBound+= w*cyvar[KEPSols2ndStage.back().cycles[i]];
    }
    
    for (int i = 0; i < KEPSols2ndStage.back().chains.size(); i++){
        int w = Chains2ndStage[Chains2ndTo3rd[KEPSols2ndStage.back().chains[i]]].AccumWeight;
        exprBound+= w*chvar[KEPSols2ndStage.back().chains[i]];
    }
    
    //cout << Beta.getName() << ">=" << exprBound << endl;
    string name = "Const11b."  + to_string(vBoundConstraint.getSize());
    const char* cName = name.c_str();
   
    vBoundConstraint.add(IloRange(env, 0, Beta - exprBound , IloInfinity, cName));
    
}
void Problem::Const11c(vector<KEPSol>&KEPSols2ndS){
    
    for (int i = 0; i < KEPSols2ndS.back().cycles.size(); i++){
        exprVxtArcsCY = IloExpr (env, 0);
        int u = KEPSols2ndS.back().cycles[i];
        for (int j = 0; j < Cycles2ndStage[Cycles2ndTo3rd[u]].get_c().size(); j++){
            exprVxtArcsCY += vertex[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]];
            if (j <= Cycles2ndStage[Cycles2ndTo3rd[u]].get_c().size() - 2){
                int v = mapArcs[make_pair(Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j], Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j + 1])];
                exprVxtArcsCY += arc[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]][v];
            }
            else{
                int v = mapArcs[make_pair(Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j], Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[0])];
                exprVxtArcsCY += arc[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]][v];
            }
        }
        //cout << cyvar[u] + exprVxtArcsCY << ">=" << 1 << endl;
        string name = "Const11c.CY"  + to_string(vFailedMatches.getSize());
        const char* cName = name.c_str();
        vFailedMatches.add(IloRange(env, 1, cyvar[u] + exprVxtArcsCY , IloInfinity, cName));
        GrandSubProb.add(vFailedMatches[vFailedMatches.getSize() - 1]);
        exprVxtArcsCY.end();
    }
    
    for (int i = 0; i < KEPSols2ndS.back().chains.size(); i++){
        exprVxtArcsCH = IloExpr (env, 0);
        int u = KEPSols2ndS.back().chains[i];
        for (int j = 0; j < Chains2ndStage[Chains2ndTo3rd[u]].Vnodes.size(); j++){
            exprVxtArcsCH += vertex[Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex];
            if (j <= Chains2ndStage[Chains2ndTo3rd[u]].Vnodes.size() - 2){
                int v = mapArcs[make_pair(Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex, Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j + 1].vertex)];
                exprVxtArcsCH += arc[Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex][v];
            }
        }
        //cout << chvar[u] + exprVxtArcsCH << ">=" << 1 << endl;
        string name = "Const11c.CH"  + to_string(vFailedMatches.getSize());
        const char* cName = name.c_str();
        vFailedMatches.add(IloRange(env, 1, chvar[u] + exprVxtArcsCH , IloInfinity, cName));
        GrandSubProb.add(vFailedMatches[vFailedMatches.getSize() - 1]);
        exprVxtArcsCH.end();
    }
    
}
bool sortChains2(Chain& c1, Chain& c2){
    return (c1.AccumWeight > c2.AccumWeight);
}
bool sortCycles2(Cycles& c1, Cycles& c2){
    return (c1.get_Many() > c2.get_Many());
}
bool CheckRepetition(vector<int>& v, map<int,bool>& included){
    for (int i = 0; i < v.size(); i++){
        map<int,bool>::iterator it;
        it = included.find(v[i]);
        if (it != included.end()){
            return true;
        }
    }
    return false;
}
bool CheckRepetition2(vector<vChain>& v, map<int,bool>& included){
    for (int i = 0; i < v.size(); i++){
        map<int,bool>::iterator it;
        it = included.find(v[i].vertex);
        if (it != included.end()){
            return true;
        }
    }
    return false;
}
bool IsSameCycle2(vector<int> c1, vector<int> c2){
    if (c1.size() == c2.size()){
        for (int j = 0; j < c1.size(); j++){
            if (c1[j] != c2[j]) return false;
        }
    }
    else{
        return false;
    }
    return true;
}
bool IsSameChain2(vector<int> c1, vector<vChain> c2){
    if (c1.size() == c2.size()){
        for (int j = 0; j < c1.size(); j++){
            if (c1[j] != c2[j].vertex) return false;
        }
    }
    else{
        return false;
    }
    return true;
}
void Problem::SampleCols2ndStage2(vector<Chain>& Chains, vector<Cycles>&Cycles, vector<IndexGrandSubSol>&SolFirstStage){
    //Organize vector by AccumWeight and HowMany
    sort(Chains.begin(), Chains.end(), sortChains2);
    sort(Cycles.begin(), Cycles.end(), sortCycles2);
    ChainNodeTPH.clear();
    CycleNodeTPH.clear();
    NodeCyChTPH.clear();
    ArcsinChCyTHP.clear();
    ArcsinChainsTHP.clear();
    ArcsinCyclesTHP.clear();
    
    //Fill in ChainNodeTPH
    for (int i = 0; i < Chains.size(); i++){
        for (int j = 0; j < Chains[i].Vnodes.size(); j++){
            ChainNodeTPH[Chains[i].Vnodes[j].vertex].push_back(i);
            NodeCyChTPH[Chains[i].Vnodes[j].vertex].push_back(i);
            if (j <= Chains[i].Vnodes.size() - 2){
                ArcsinChainsTHP[make_pair(Chains[i].Vnodes[j].vertex, Chains[i].Vnodes[j + 1].vertex)].push_back(i);
                ArcsinChCyTHP[make_pair(Chains[i].Vnodes[j].vertex, Chains[i].Vnodes[j + 1].vertex)].push_back(i);
            }
        }
    }

    //Fill in CycleNodeTPH
    for (int i = 0; i < Cycles.size(); i++){
        for (int j = 0; j < Cycles[i].get_c().size(); j++){
            CycleNodeTPH[Cycles[i].get_c()[j]].push_back(i);
            NodeCyChTPH[Cycles[i].get_c()[j]].push_back(i);
            if (j <= Cycles[i].get_c().size() - 2){
                ArcsinCyclesTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[j + 1])].push_back(i);
                ArcsinChCyTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[j + 1])].push_back(i);
            }
            else{
                ArcsinCyclesTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[0])].push_back(i);
                ArcsinChCyTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[0])].push_back(i);
            }
        }
    }
    
    CycleNodeSPH.clear();
    ChainNodeSPH.clear();
    
    //Find vertex-disjoint solution
    int counter = 0;
    map<int,bool>included;
    KEPSols2ndStage.push_back(KEPSol());
    
    for (int j = 0; j < SolFirstStage.size(); j++){
        if (SolFirstStage[j].get_cc()[0] < Pairs){// It is a cycle
            for (int i = Cycles.size() - 1; i >= 0; i--){
                //If Cycle[i] is in SolFirstStage
                if (IsSameCycle2(SolFirstStage[j].get_cc(), Cycles[i].get_c()) == true){
                    Cycles2ndTo3rd[counter] = i;
                    Cycles3rdTo2nd[i] = counter;
                    KEPSols2ndStage.back().cycles.push_back(counter);
                    counter++;
                    break;//cycle found
                }
            }
        }
    }

    counter = 0;
    for (int j = 0; j < SolFirstStage.size(); j++){
        if (SolFirstStage[j].get_cc()[0] >= Pairs){// It is a chain
            for (int i = 0; i < Chains.size(); i++){
                //If Cycle[i] is in SolFirstStage
                if (IsSameChain2(SolFirstStage[j].get_cc(), Chains[i].Vnodes) == true){
                    Chains2ndTo3rd[counter] = i;
                    Chains3rdTo2nd[i] = counter;
                    KEPSols2ndStage.back().chains.push_back(counter);
                    counter++;
                    break;
                }
            }
        }
    }
    
    cout << endl;
}
