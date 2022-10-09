//
//  GrandSubProbMasterPICEF.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2022-05-05.
//  Copyright © 2022 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandSubProbMasterPICEF.hpp"


void Problem::GrandSubProbMasterPICEF(vector<Cycles>&Cycles2ndStage, vector<IndexGrandSubSol>&SolFirstStage){
    // Create model
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    GrandSubProb = IloModel(env);
    cplexGrandSubP = IloCplex(GrandSubProb);
    cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexGrandSubP.setParam(IloCplex::Param::Threads, 1);
    cplexGrandSubP.setOut(env.getNullStream());
    tStart2ndS = clock();
    
    SampleCols2ndStagePICEF(Cycles2ndStage, SolFirstStage);
    
    //Get selected vertices
    vector<int>ListSelVertices = GetSelPairs(SolFirstStage);
    
    //Create cycle variables
    cyvar = IloNumVarArray(env, Cycles2ndTo3rd.size(), 0, 1, ILOINT);
    for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
       SetName1Index(cyvar[i], "x", i + 1);
       //cout << r[i].getName() << endl;
    }
    
    //Create chain variables
    E_sijl = Create4DBinTHP(env, Chains3rdStageP.back(), int(ChainLength), "E_sijl");
    
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
    Const13b(KEPSols2ndStage, Chains3rdStageP.back(), ListSelVertices);
    
    //Create Active Cols Constraint
    vFailedMatches = IloRangeArray(env);
    Const13c(KEPSols2ndStage);
    
    vFailedArcZero = IloRangeArray(env);
    Const13d(Chains3rdStageP.back(), E_sijl);
    
    vFailedArcs = IloRangeArray(env);
    Const13e(Chains3rdStageP.back(), E_sijl);
    
    //Do not pick a failed arc adjacent to a failed vertex
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        IloExpr in (env,0);  IloExpr out (env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            out+= arc[i][j];
        }
        for (int j = 0; j < PredMap[i].size(); j++){
            in+= arc[PredMap[i][j].first][PredMap[i][j].second];
        }
        GrandSubProb.add(IloRange(env, -IloInfinity, in + out + 2*vertex[i], 2));
        in.end();
        out.end();
    }
    
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
    
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    cplexGrandSubP.solve();
    tTotalMP2ndPH+= (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
    if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
        cout << "S.O.S. This should not happen." << endl;
    }
    else{
        //cout << cplexGrandSubP.getStatus();
        SPMIP_Obj = cplexGrandSubP.getValue(Beta);
        vertex_sol = IloNumArray(env, Nodes);
        cplexGrandSubP.getValues(vertex_sol,vertex);
        
        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            cplexGrandSubP.getValues(arc_sol[f],arc[f]);
        }
        
        //Solve 2nd Stage
        BendersPICEF(Cycles2ndStage, SolFirstStage);
        tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
        GlobalIte2ndStage += Ite2ndS;
        //Send scenario to first phase
        if (Ite1stStage == 1){
            scenarios.clear();
        }
        scenarios.push_back(map<pair<int,int>,bool>());
        for (auto it = OptFailedArcs.begin(); it != OptFailedArcs.end(); it++){
            scenarios.back()[it->first] = true;
        }
        for (auto it = OptFailedVertices.begin(); it != OptFailedVertices.end(); it++){
            scenarios.back()[make_pair(-1, it->first)] = true;
        }

        tStart1stS = clock();
        if (Ite1stStage == 1){
            //Update Con7g
            for (int j = 0; j < Pairs; j++) { // para todo j
                auto it = scenarios[0].find(make_pair(-1, j));
                if (it != scenarios[0].end()){
                    for (int c = 0; c < CycleNode[j].size(); c++) {
                        X_cu[CycleNode[j][c]][0].setUB(0);
                    }
                    for (int l = 0; l < ChainLength; l++) {
                        for (int i = 0; i < PredMap[j].size(); i++) {
                            int i2 = PredMap[j][i].first;
                            int j2 = PredMap[j][i].second;
                            //cout << E_ijlu[i2][j2][l][Ulast].getName();
                            E_ijlu[i2][j2][l][0].setUB(0); // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
                        }
                    }
                }
            }
            //Update Con7h
            for (int j = Pairs; j < AdjacencyList.getSize(); j++) { // para todo j
                auto it = scenarios[0].find(make_pair(-1, j));
                if (it != scenarios[0].end()){
                    for (int i = 0; i < AdjacencyList[j].getSize(); i++) {
                        E_ijlu[j][i][0][0].setUB(0);
                    }
                }
            }
            //Constraint (7j)
            if (MaxArcFailures > 0){
                //Constraint (7k)
                CreateCon7k(env, X_cu, E_ijlu, scenarios, 0, AdjacencyList, CycleArcs, PredMap, int(ChainLength), "7k");
            }
        }
        else{
            /////////////////////Create new cols/////////////////////////
            // Y_ju
            for (int i = 0; i < Y_ju.getSize(); i++){
                //Create new Y_ju column
                IloNumColumn col(env);
                //Create new Y_ju variable
                Y_ju[i].add(IloIntVar(col, 0, 1));
                string varName = "Y_ju_" + to_string(i) + "_" + to_string(scenarios.size() - 1);
                Y_ju[i][Y_ju[i].getSize() - 1].setName(varName.c_str());
            }
            //X_cu
            for (int i = 0; i < X_cu.getSize(); i++){
                //Create new X_cu column
                IloNumColumn col(env);
                //Create new Y_ju variable
                X_cu[i].add(IloIntVar(col, 0, 1));
                string varName = "X_cu_" + to_string(i) + "_" + to_string(scenarios.size() - 1);
                X_cu[i][X_cu[i].getSize() - 1].setName(varName.c_str());
            }
            // E_ijlu
            for (int i = 0; i < E_ijlu.getSize(); i++){
                for (int j = 0; j < E_ijlu[i].getSize(); j++){
                    for (int l = 0; l < E_ijlu[i][j].getSize(); l++){
                        //Create new E_ijlu column
                        IloNumColumn col(env);
                        //Create new Y_ju variable
                        E_ijlu[i][j][l].add(IloIntVar(col, 0, 1));
                        string varName = "E_ijlu_" + to_string(i) + "_" + to_string(int(AdjacencyList[i][j] - 1)) + "_" + to_string(l) + "_" + to_string(scenarios.size() - 1);
                        E_ijlu[i][j][l][scenarios.size() - 1].setName(varName.c_str());
                        //cout << E_ijlu[i][j][l][scenarios.size() - 1].getName() << endl;
                    }
                }
            }
            SetUB4DBin (env, E_ijlu, distNDD, int(Pairs), int(ChainLength), int (scenarios.size() - 1));
            
            /////////////////////Call constraints/////////////////////////
            //Constraint (7b)
            IloRangeArray cons7b(env);
            cons7b = CreateCons7b(env, Z, Y_ju, int(scenarios.size() - 1), int(Pairs), "7b");
            RobustMod.add(cons7b);
            //Constraint (7c)
            IloRangeArray cons7c (env);
            cons7c = CreateCons7c(env, Y_ju, X_c, E_ijl, int(scenarios.size() - 1), int(Pairs), CycleNode, PredMap, "7c");
            RobustMod.add(cons7c);
            //Constraint (7d)
            IloRangeArray cons7d (env);
            cons7d = CreateCon7d(env, Y_ju, X_cu, E_ijlu, int(scenarios.size() - 1), int(Pairs), CycleNode, int(ChainLength), PredMap, "7d");
            RobustMod.add(cons7d);
            //Constraint (7g)
            IloRangeArray cons7g (env);
            cons7g = CreateCon7g(env, X_cu, E_ijlu, scenarios, int(scenarios.size() - 1), int(Pairs), CycleNode, PredMap, int(ChainLength), "7g");
            RobustMod.add(cons7g);
            //Constraint (7h)
            IloRangeArray cons7h (env);
            cons7h = CreateCon7h(env, E_ijlu, scenarios, int(scenarios.size() - 1), int(Pairs), AdjacencyList, "7h");//Set UB to 0 when U_uj = 1
            RobustMod.add(cons7h);
            //Constraint (7j)
            IloRangeArray cons7j (env);
            cons7j = CreateCon7j(env, E_ijlu, int(scenarios.size() - 1), int(Pairs), int(ChainLength), PredMap, AdjacencyList, "7j");
            RobustMod.add(cons7j);
            if (MaxArcFailures > 0){
                //Constraint (7k)
                CreateCon7k(env, X_cu, E_ijlu, scenarios, int(scenarios.size() - 1), AdjacencyList, CycleArcs, PredMap, int(ChainLength), "7k");
            }
        }
        
        //cplexRobust.exportModel("RO_Model.lp");
        LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
        if (LeftTime < 0){
            GlobalIte2ndStage += Ite2ndS;
            tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            Print2ndStage("TimeOut", SolFirstStage);
        }
        cplexRobust.setParam(IloCplex::Param::TimeLimit, LeftTime);
        cplexRobust.solve();
        tTotal1stS += (clock() - tStart1stS)/double(CLOCKS_PER_SEC);
        if (cplexRobust.getStatus() == IloAlgorithm::Infeasible){
            cout << "This should not happen." << endl;
        }
        else if (cplexRobust.getStatus() == IloAlgorithm::Optimal){
            //////////////////Retrieve solution/////////////////////////
            FPMIP_Obj = cplexRobust.getObjValue();
            vector<IndexGrandSubSol>SolFirstStageNew;
            IloNumArray xsol(env, ListCycles.size());
            cplexRobust.getValues(xsol,X_c);
            for (int i = 0; i < xsol.getSize(); i++){
                if (xsol[i] > 0.9){
                    SolFirstStageNew.push_back(IndexGrandSubSol(ListCycles[i].get_c(), ListCycles[i].get_c().size()));
                }
            }
            IloNumArray3 esol(env, AdjacencyList.getSize());
            for (int i = 0; i < esol.getSize(); i++){
                esol[i] = IloNumArray2 (env, AdjacencyList[i].getSize());
                for (int j = 0; j < esol[i].getSize(); j++){
                    esol[i][j] = IloNumArray(env, ChainLength);
                    cplexRobust.getValues(esol[i][j],E_ijl[i][j]);
                    for (int k = 0; k < E_ijl[i][j].getSize(); k++){
                        //if (esol[i][j][k] > 0.9) cout << E_ijl[i][j][k].getName() << "\t";
                    }
                }
            }
            vector<vector<int>>vChains;
            vChains = GetChainsFrom1stStageSol(AdjacencyList,esol, Pairs, ChainLength);
            for (int i = 0; i < vChains.size(); i++){
                SolFirstStageNew.push_back(IndexGrandSubSol(vChains[i], vChains[i].size() - 1));
            }
            if (SPMIP_Obj < FPMIP_Obj){ //Send scenario to 1st. stage
                //Call 2nd. stage
                tStart2ndS = clock();
                Cycles2ndStage = Get2ndStageCycles (SolFirstStageNew, RecoursePolicy);
                tTotalFindingCyCh+= (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                GrandSubProb.end();
                cplexGrandSubP.end();
                ub_tcyvar.clear();
                SPMIP_Obj = 0;
                mTHPMIP.end();
                cplexmTHPMIP.end();
                Ite2ndS = 0;
                Ite1stStage++;
                LOWEST_TPMIP_Obj = INT_MAX;
                GrandSubProbMasterPICEF(Cycles2ndStage, SolFirstStage);
            }
            else if (SPMIP_Obj == FPMIP_Obj){
                //Robust Solution found
                Print2ndStage("Optimal", SolFirstStage);
            }
            else{
                cout << endl << "This should never happen";
            }
        }
        else{
            GlobalIte2ndStage += Ite2ndS;
            tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            Print2ndStage("TimeOut", SolFirstStage);
        }
        
    }

}
void Problem::Const13b(vector<KEPSol>&KEPSols2ndStage, vector<IndexGrandSubSol>&ChainsTPH, vector<int>&SelectedVertices){
    
    exprBound = IloExpr(env,0);
    for (int i = 0; i < KEPSols2ndStage.back().cycles.size(); i++){
        int w = Cycles2ndStage[Cycles2ndTo3rd[KEPSols2ndStage.back().cycles[i]]].get_Many();
        exprBound+= w*cyvar[KEPSols2ndStage.back().cycles[i]];
    }
    
    for (int i = 0; i < ChainsTPH.size(); i++){
        for (int j = 0; j < ChainsTPH[i].get_cc().size() - 1; j++){
            int w = 0;
            if (wasSelected(SelectedVertices, ChainsTPH[i].get_cc()[j + 1]) == true){
                w = 1;
            }
            exprBound+= w*E_sijl[E_sijl.getSize() - 1][i][j][j];
        }
    }
    
    //cout << Beta.getName() << ">=" << exprBound << endl;
    string name = "Const13b."  + to_string(vBoundConstraint.getSize());
    const char* cName = name.c_str();
   
    vBoundConstraint.add(IloRange(env, 0, Beta - exprBound , IloInfinity, cName));
    GrandSubProb.add(vBoundConstraint[vBoundConstraint.getSize() - 1]);
    exprBound.end();
        
}
void Problem::Const13c(vector<KEPSol>&KEPSols2ndS){
    
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
        string name = "Const13c.CY"  + to_string(vFailedMatches.getSize());
        const char* cName = name.c_str();
        vFailedMatches.add(IloRange(env, 1, cyvar[u] + exprVxtArcsCY , IloInfinity, cName));
        GrandSubProb.add(vFailedMatches[vFailedMatches.getSize() - 1]);
        exprVxtArcsCY.end();
    }
    
}
void Problem::Const13d(vector<IndexGrandSubSol>&ChainsTPH, NumVar4D& E_sijl){
    for (int s = 0; s < ChainsTPH.size(); s++){
        int v = mapArcs[make_pair(ChainsTPH[s].get_cc()[0], ChainsTPH[s].get_cc()[1])];
        string name = "Const13d"  + to_string(vFailedArcZero.getSize());
        const char* cName = name.c_str();
        
        //cout << E_sijl[E_sijl.getSize() - 1][s][0][0].getName() << '\t' << vertex[ChainsTPH[s].get_cc()[0]].getName() << '\t' << vertex[ChainsTPH[s].get_cc()[1]].getName() << '\t' << arc[ChainsTPH[s].get_cc()[0]][v].getName() << endl;
        
        vFailedArcZero.add(IloRange(env, 1, E_sijl[E_sijl.getSize() - 1][s][0][0] + vertex[ChainsTPH[s].get_cc()[0]] + vertex[ChainsTPH[s].get_cc()[1]] + arc[ChainsTPH[s].get_cc()[0]][v], IloInfinity, cName));
        GrandSubProb.add(vFailedArcZero[vFailedArcZero.getSize() - 1]);
    }
}
void Problem::Const13e(vector<IndexGrandSubSol>& ChainsTPH, NumVar4D& E_sijl){
    for (int s = 0; s < ChainsTPH.size(); s++){
        for (int i = 1; i < ChainsTPH[s].get_cc().size() - 1; i++){
            int v = mapArcs[make_pair(ChainsTPH[s].get_cc()[i], ChainsTPH[s].get_cc()[i + 1])];
            string name = "Const13e"  + to_string(vFailedArcs.getSize());
            const char* cName = name.c_str();
            
            //cout << E_sijl[E_sijl.getSize() - 1][s][i][i].getName() << '\t' << E_sijl[E_sijl.getSize() - 1][s][i - 1][i - 1].getName() << '\t' << vertex[ChainsTPH[s].get_cc()[i + 1]].getName() << '\t' << arc[ChainsTPH[s].get_cc()[i]][v].getName() << endl;

            
            vFailedArcs.add(IloRange(env, 0, E_sijl[E_sijl.getSize() - 1][s][i][i] - E_sijl[E_sijl.getSize() - 1][s][i - 1][i - 1] + vertex[ChainsTPH[s].get_cc()[i + 1]] +  arc[ChainsTPH[s].get_cc()[i]][v], IloInfinity, cName));
            GrandSubProb.add(vFailedArcs[vFailedArcs.getSize() - 1]);
        }
    }
}
bool wasSelected(vector<int>&v, int a){
    for (int i = 0; i < v.size(); i++) {
        if (v[i] == a) return true;
    }
    return false;
}
void Problem::SampleCols2ndStagePICEF(vector<Cycles>&Cycles, vector<IndexGrandSubSol>&SolFirstStage){
    
    //Find vertex-disjoint solution
    int counter = 0;
    map<int,bool>included;
    KEPSols2ndStage.clear();
    Chains3rdStageP.clear();
    Cycles2ndTo3rd.clear();
    Cycles3rdTo2nd.clear();
    CycleNodeTPH.clear();
    KEPSols2ndStage.push_back(KEPSol());
    
    //Fill in CycleNodeTPH
    for (int i = 0; i < Cycles.size(); i++){
        for (int j = 0; j < Cycles[i].get_c().size(); j++){
            CycleNodeTPH[Cycles[i].get_c()[j]].push_back(i);
        }
    }
    
    
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
    
    Chains3rdStageP.push_back(vector<IndexGrandSubSol>());
    for (int j = 0; j < SolFirstStage.size(); j++){
        if (SolFirstStage[j].get_cc()[0] >= Pairs){// It is a chain
            Chains3rdStageP.back().push_back(SolFirstStage[j]);
        }
    }
    
}
NumVar4D Create4DBinTHP(IloEnv& env, vector<IndexGrandSubSol>& ChainsTHP, int L, string name){
    NumVar4D chvar(env, 1);
    chvar[0] = NumVar3D(env, ChainsTHP.size());
    for (int s = 0; s < ChainsTHP.size(); s++){
        chvar[0][s] = NumVar2D(env, ChainsTHP[s].get_cc().size() - 1);
        for (int i = 0; i < ChainsTHP[s].get_cc().size() - 1; i++){
            chvar[0][s][i] = IloNumVarArray(env, L, 0, 1, ILOINT);
            for (int j = 0; j < L; j++){
                if (i != j){
                    chvar[0][s][i][j].setUB(0);
                }
                else{
                    string auxName = name + "_" + to_string(chvar.getSize() - 1) + "_" + to_string(int(ChainsTHP[s].get_cc()[i] + 1)) + "_" + to_string(int(ChainsTHP[s].get_cc()[i + 1] + 1)) + "_" + to_string(int(j));
                    chvar[0][s][i][j].setName(auxName.c_str());
                }
            }
        }
    }
    
    return chvar;
}
IloRangeArray CreateCon18d(IloEnv& env, NumVar2D& var, IloNumArray2& AdjacencyList, map<int, vector<pair<int,int>>>& PredMap, int Pairs){
    IloRangeArray range (env);
    
    for (int j = 0 ; j < AdjacencyList.getSize(); j++){
        if (j < Pairs){
            IloExpr expr (env, 0);
            expr += IloSum(var[j]);
            for (int i = 0; i < PredMap[j].size(); i++) {
                int i2 = PredMap[j][i].first;
                int j2 = PredMap[j][i].second;
                expr-= var[i2][j2];
            }
            //cout << endl << expr << endl;
            string name2 = "const18d[" + to_string(j) + "]";
            range.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
        }
    }
    
    return range;
}
IloRangeArray CreateCon18f(IloEnv& env, NumVar2D& var2D, NumVar3D& var3D){
    IloRangeArray range (env);
    
    for (int i = 0 ; i < var3D.getSize(); i++){
        for (int j = 0; j < var3D[i].getSize(); j++){
            string name2 = "const18f[" + to_string(i) + "][" + to_string(j) + "]";
            //cout << var2D[i][j].getName() << "\t" << IloSum(var3D[i][j]) << endl;
            range.add(IloRange(env, -IloInfinity, var2D[i][j] - IloSum(var3D[i][j]), 0, name2.c_str()));
        }
    }
    
    return range;
}
void GetUB_yij(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices, NumVar2D& var2D, IloNumArray2& AdjaList){
    
    for (int i = 0; i < var2D.getSize(); i++){
        for (int j = 0; j < var2D[i].getSize(); j++){
            auto it = FailedVertices.find(i);
            if (it == FailedVertices.end()){
                it = FailedVertices.find(AdjaList[i][j] - 1);
                if (it == FailedVertices.end()){
                    auto it2 = FailedArcs.find(make_pair(i, AdjaList[i][j] - 1));
                    if (it2 == FailedArcs.end()){
                        var2D[i][j].setUB(1);
                    }
                    else{
                        var2D[i][j].setUB(0);
                        //cout << var2D[i][j].getName() << endl;
                    }
                }
                else{
                    var2D[i][j].setUB(0);
                    //cout << var2D[i][j].getName() << endl;
                }
            }
            else{
                var2D[i][j].setUB(0);
                //cout << var2D[i][j].getName() << endl;
            }
        }
    }

}
void Problem::BendersPICEF(vector<Cycles>& Cycles2ndStage, vector<IndexGrandSubSol>&SolFirstStage){
    tStartRecoMIP = clock();
    //Get scenario
    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
    //Get selected vertices
    vector<int>SelectedVertices = GetSelPairs(SolFirstStage);
    //Create model
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    if (LeftTime < 0){
        GlobalIte2ndStage += Ite2ndS;
        tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
        Print2ndStage("TimeOut", SolFirstStage);
    }
    mTHPMIP = IloModel (env);
    cplexmTHPMIP = IloCplex(mTHPMIP);
    cplexmTHPMIP.setParam(IloCplex::Param::Threads, 1);
    cplexmTHPMIP.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexmTHPMIP.setOut(env.getNullStream());
    
    //Create decision variables and set decision variables values according to the scenario
    tcyvar = Create_tcyvar("r", Cycles2ndStage);
    
    yij = NumVar2D(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        yij[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, 1, ILOINT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            SetName2(yij[i][j], "y", i + 1, AdjacencyList[i][j]);
        }
    }
    
    Eijl = Create3DBin (env, AdjacencyList, int(ChainLength), "Eijl");
    
    SetUB3DBin (env, Eijl, distNDD, int(Pairs), int(ChainLength));
    
    //Add objective value
    //Change variables' UB in 3rd phase
    ub_tcyvar.clear(); // <cycle number, 0 or 1>
    ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
    
    IloExpr obj(env,0);
    obj = GetObjTPH_PICEF(env, Cycles2ndStage, SelectedVertices, tcyvar, yij, Eijl, AdjacencyList, int(Pairs), int(ChainLength), ub_tcyvar);
    string name = "Obj_THP_PICEF";
    ObjTHP = IloObjective(env, obj, IloObjective::Maximize, name.c_str());
    mTHPMIP.add(ObjTHP);
    
    // Constraint (18b)
    IloRangeArray cons18b (env);
    cons18b = CreateCon7e(env, tcyvar, Eijl, int(Pairs), CycleNodeTPH, PredMap, "cons18b");
    mTHPMIP.add(cons18b);
    
    //Constraint (18c)
    IloRangeArray cons18c (env);
    cons18c = CreateCon7f(env, Eijl, int(Pairs), AdjacencyList, "cons18c");
    mTHPMIP.add(cons18c);
    
    //Constraint (18d)
    IloRangeArray cons18d (env);
    cons18d = CreateCon18d(env, yij, AdjacencyList, PredMap, int(Pairs));
    mTHPMIP.add(cons18d);
    
    //Constraint (18f)
    IloRangeArray cons18f (env);
    cons18f = CreateCon18f(env, yij, Eijl);
    mTHPMIP.add(cons18f);
    
    //Constraint (18g)
    IloRangeArray cons18g (env);
    cons18g = CreateCon7i(env, Eijl, Pairs, int(ChainLength), PredMap, AdjacencyList, "cons18g");
    mTHPMIP.add(cons18g);
    
    while(true){
        Ite2ndS++;
        GetUB_yij(FailedArcs, FailedVertices, yij, AdjacencyList);
        //cplexmTHPMIP.exportModel("THPMIP.lp");
        LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
        if (LeftTime < 0){
            GlobalIte2ndStage += Ite2ndS;
            tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            Print2ndStage("TimeOut", SolFirstStage);
        }
        cplexmTHPMIP.setParam(IloCplex::Param::TimeLimit, LeftTime);
        cplexmTHPMIP.solve();
        tTotalRecoMIP += (clock() - tStartRecoMIP)/double(CLOCKS_PER_SEC);
        if (cplexmTHPMIP.getStatus() == IloAlgorithm::Infeasible){
            cout << "S.O.S. This should not happen." << endl;
        }
        else if (cplexmTHPMIP.getStatus() == IloAlgorithm::Optimal){
            tStartRecoMIP = clock();
            TPMIP_Obj = cplexmTHPMIP.getObjValue();
            double Actual_TPHObj = 0;
            IloNumArray xsol(env, Cycles2ndStage.size());
            cplexmTHPMIP.getValues(xsol,tcyvar);
            KEPSols2ndStage.push_back(KEPSol());
            vector<KEPSol> KEPSol2ndStageMissing;
            KEPSol2ndStageMissing.push_back(KEPSol());
            
            for (int i = 0; i < xsol.getSize(); i++){
                if (xsol[i] > 0.9){
                    auto find = ub_tcyvar.find(i);
                    if (find == ub_tcyvar.end())    Actual_TPHObj += Cycles2ndStage[i].get_Many();
                    auto it = Cycles3rdTo2nd.find(i);
                    if (it == Cycles3rdTo2nd.end()){
                        int n = int(Cycles2ndTo3rd.size());
                        Cycles2ndTo3rd[n] = i;
                        Cycles3rdTo2nd[i] = n;
                        KEPSols2ndStage.back().cycles.push_back(Cycles3rdTo2nd[i]);
                        KEPSol2ndStageMissing.back().cycles.push_back(Cycles3rdTo2nd[i]);
                    }
                    else{
                        KEPSols2ndStage.back().cycles.push_back(Cycles3rdTo2nd[i]);
                    }
                }
            }
            
            IloNumArray2 ysol(env, AdjacencyList.getSize());
            // cout << endl << "sol" << endl;
            for (int i = 0; i < ysol.getSize(); i++){
                ysol[i] = IloNumArray (env, AdjacencyList[i].getSize());
                cplexmTHPMIP.getValues(ysol[i], yij[i]);
                //for (int j = 0; j < ysol[i].getSize(); j++){
                //    if (ysol[i][j] > 0.9){
                //        cout << yij[i][j].getName() << "\t";
                //    }
                //}
            }
        //    cout << endl << "Chains: " << endl;
            IloNumArray3 esol(env, AdjacencyList.getSize());
            for (int i = 0; i < esol.getSize(); i++){
                esol[i] = IloNumArray2 (env, AdjacencyList[i].getSize());
                for (int j = 0; j < esol[i].getSize(); j++){
                    esol[i][j] = IloNumArray(env, ChainLength);
                    cplexmTHPMIP.getValues(esol[i][j],Eijl[i][j]);
                    //for (int l = 0; l < ChainLength; l++){
                    //    if (esol[i][j][l] > 0.9){
                    //        cout << Eijl[i][j][l].getName() << "\t";
                    //    }
                    //}
                }
            }
            
            vector<vector<int>>vChains;
            vChains = GetChainsFrom1stStageSol(AdjacencyList,esol, int(Pairs), int(ChainLength));
            
            Chains3rdStageP.push_back(vector<IndexGrandSubSol>());
            for (int i = 0; i < vChains.size(); i++){
                for (int j = 0; j < vChains[i].size() - 1; j++){
                    int v = mapArcs[make_pair(vChains[i][j], vChains[i][j + 1])];
                    if (yij[vChains[i][j]][v].getUB() != 0){
                        //cout << yij[vChains[i][j]][v].getUB() << " " << yij[vChains[i][j]][v].getName() << endl;
                        if (wasSelected(SelectedVertices, vChains[i][j + 1]) == true){
                            Actual_TPHObj++;
                        }
                    }
                    else{
                        break;
                    }
                }
                Chains3rdStageP.back().push_back(IndexGrandSubSol(vChains[i], 0));
            }
            
            TPMIP_Obj = Actual_TPHObj;
            tTotalRecoMIP += (clock() - tStartRecoMIP)/double(CLOCKS_PER_SEC);
            tStartMP2ndPH = clock();
            
            
            if (TPMIP_Obj - SPMIP_Obj <= 0.0002){
                OptFailedArcs = FailedArcs;
                OptFailedVertices = FailedVertices;
                break;
            }
            else if (SPMIP_Obj <= TPMIP_Obj){
                //Create new cycle columns 2ndPhase
                for (int i = 0; i < KEPSol2ndStageMissing.back().cycles.size(); i++){
                    int c = KEPSol2ndStageMissing.back().cycles[i];
                    //Create new cycle column
                    IloNumColumn col(env);
                    //Create new cycle variable
                    cyvar.add(IloIntVar(col, 0, 1));
                    string name = "x." + to_string(c);
                    const char* varName = name.c_str();
                    cyvar[cyvar.getSize() - 1].setName(varName);
                    cyvar[cyvar.getSize() - 1].setBounds(0, 1);
                }
                //Create new chain columns 2ndPhase
                E_sijl[E_sijl.getSize() - 1].add(NumVar3D(env, Chains3rdStageP.back().size()));
                for (int i = 0; i < Chains3rdStageP.back().size(); i++){
                    E_sijl[E_sijl.getSize() - 1][i] = NumVar2D(env, Chains3rdStageP.back()[i].get_cc().size() - 1);
                    for (int j = 0; j < Chains3rdStageP.back()[i].get_cc().size() - 1; j++){
                        E_sijl[E_sijl.getSize() - 1][i][j] = IloNumVarArray(env);
                        for (int l = 0; l < ChainLength; l++){
                            //Create new cycle column
                            IloNumColumn col(env);
                            //Create new chain variable
                            E_sijl[E_sijl.getSize() - 1][i][j].add(IloIntVar(col, 0, 1));
                            string auxName = "E_sijl_" + to_string(E_sijl.getSize() - 1) + "_" + to_string(int(Chains3rdStageP.back()[i].get_cc()[j] + 1)) + "_" + to_string(int(Chains3rdStageP.back()[i].get_cc()[j + 1] + 1)) + "_" + to_string(int(j));
                            E_sijl[E_sijl.getSize() - 1][i][j][l].setName(auxName.c_str());
                            if (j != l) E_sijl[E_sijl.getSize() - 1][i][j][l].setUB(0);
                        }
                    }
                }
                
                //Create constraints 2ndPhase
                //Constraint (13b)
                Const13b(KEPSols2ndStage, Chains3rdStageP.back(), SelectedVertices);
                
                //Constraint (13c)
                Const13c(KEPSol2ndStageMissing);
                
                //Constraint (13d)
                Const13d(Chains3rdStageP.back(), E_sijl);
                
                //Constraint (13e)
                Const13e(Chains3rdStageP.back(), E_sijl);
                
                //Resolve 2ndPhase
                    // set left time
                //cplexGrandSubP.exportModel("GrandSubP.lp");
                LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
                if (LeftTime < 0){
                    GlobalIte2ndStage += Ite2ndS;
                    tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                    Print2ndStage("TimeOut", SolFirstStage);
                }
                cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, LeftTime);
                cplexGrandSubP.solve();
                tTotalMP2ndPH += (clock() - tStartMP2ndPH)/double(CLOCKS_PER_SEC);
                if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
                    cout << "S.O.S. This should not happen." << endl;
                }
                else if (cplexGrandSubP.getStatus() == IloAlgorithm::Optimal){
                    tStartMP2ndPH = clock();
                    //cout << cplexGrandSubP.getStatus();
                    SPMIP_Obj = cplexGrandSubP.getValue(Beta);
                    vertex_sol = IloNumArray(env, Nodes);
                    cplexGrandSubP.getValues(vertex_sol,vertex);
                    
                    arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
                    for (int f = 0; f < arc_sol.getSize(); f++){
                        arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
                        cplexGrandSubP.getValues(arc_sol[f],arc[f]);
                    }
                    //Update scenario
                    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
                    tTotalMP2ndPH += (clock() - tStartMP2ndPH)/double(CLOCKS_PER_SEC);
                    //Repeat
                }
                else{
                    GlobalIte2ndStage += Ite2ndS;
                    tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                    Print2ndStage("TimeOut", SolFirstStage);
                }
            }
            else{
                cout << "S.O.S";
            }
            tStartRecoMIP = clock();
            //Update TPH objective coefficients
            ub_tcyvar.clear(); // <cycle number, 0 or 1>
            ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
            for(int i = 0; i < tcyvar.getSize(); i++){
                map<int,bool>::iterator find = ub_tcyvar.find(i);
                if (find != ub_tcyvar.end()){
                    //Change ObjCoefficients
                    ObjTHP.setLinearCoef(tcyvar[i], 1);
                }
                else{
                    //Change ObjCoefficients
                    int n = int(Cycles2ndStage[i].get_Many());
                    ObjTHP.setLinearCoef(tcyvar[i], (n*Nodes + 1));
                }
            }
            mTHPMIP.add(ObjTHP);
        }
        else{
            GlobalIte2ndStage += Ite2ndS;
            tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            Print2ndStage("TimeOut", SolFirstStage);
        }
    }
    
}
IloExpr GetObjTPH_PICEF(IloEnv& env, vector<Cycles>&Cycles2ndStage, vector<int>&SelectedVertices, IloNumVarArray& var1D, NumVar2D& var2D, NumVar3D& var3D, IloNumArray2& AdjaList, int Pairs, int L, map<int,bool>& ub_tcyvar){
    IloExpr expr(env);

    for (int i = 0; i < Cycles2ndStage.size(); i++){
        auto it = ub_tcyvar.find(i);
        if (it == ub_tcyvar.end()){
            expr+= var1D[i]*((Cycles2ndStage[i].get_Many()*AdjaList.getSize()) + 1);
        }
        else{
            expr+= var1D[i];
        }
    }
    
    for (int i = 0; i < AdjaList.getSize(); i++){
        for (int j = 0; j < AdjaList[i].getSize(); j++){
            if (i < Pairs){
                if (wasSelected(SelectedVertices, AdjaList[i][j] - 1) == true){
                    expr+= AdjaList.getSize()*var2D[i][j];
                }
            }
            else{
                if (wasSelected(SelectedVertices, AdjaList[i][j] - 1) == true){
                    expr+= (AdjaList.getSize() + 1)*var2D[i][j];
                }
            }
        }
    }
    double n = AdjaList.getSize();
    double d = 1.0/n;
    for (int i = 0; i < AdjaList.getSize(); i++){
        for (int j = 0; j < AdjaList[i].getSize(); j++){
            for (int l = 0; l < L; l++){
                expr+= d*var3D[i][j][l];
            }
        }
    }
    
    return expr;
}
