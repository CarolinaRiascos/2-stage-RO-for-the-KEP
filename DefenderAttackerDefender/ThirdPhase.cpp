//
//  ThirdPhase.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "ThirdPhase.hpp"
void Problem::THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<int>&vinFirstStage){
    tStartRecoMIP = clock();
    //Get scenario
    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
    //Create model
    mTHPMIP = IloModel (env);
    
    //Create decision variables and set decision variables values according to the scenario
    tcyvar = Create_tcyvar("r", Cycles2ndStage);
    tchvar = Create_tchvar("s", Chains2ndStage);
    
    //Create constraints
    IloRangeArray Disjoint(env);
    Disjoint = DisjointTHP(tcyvar, tchvar);
    mTHPMIP.add(Disjoint);
    
    //Add objective function
    //Set Objective
    IloExpr obj(env,0);
    obj = GetObjTPH(Cycles2ndStage, Chains2ndStage, THP_Method);
    string name = "Obj_THP";
    ObjTHP = IloObjective(env, obj, IloObjective::Maximize, name.c_str());
    mTHPMIP.add(ObjTHP);
    tTotalRecoMIP += (clock() - tStartRecoMIP)/double(CLOCKS_PER_SEC);
    cplexmTHPMIP = IloCplex(mTHPMIP);
    cplexmTHPMIP.setParam(IloCplex::Param::Threads, 1);
    cplexmTHPMIP.setOut(env.getNullStream());
    Ite2ndS = 0;
    while(true){
        LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
        if (LeftTime <= 0){
            break;
        }
        Ite2ndS++;
        //Solve formulation
        tStartCG = clock();
        ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
        ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
        bool runColGen = ColumnGeneration(ub_tcyvar, ub_tchvar);
        tTotalRecoCG += (clock() - tStartCG)/double(CLOCKS_PER_SEC);
        if (runColGen == false){
            tStartRecoMIP = clock();
            LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
            if (LeftTime < 0) break;
            cplexmTHPMIP.setParam(IloCplex::Param::TimeLimit, LeftTime);
            cplexmTHPMIP.solve();
            tTotalRecoMIP += (clock() - tStartRecoMIP)/double(CLOCKS_PER_SEC);
            if (cplexmTHPMIP.getStatus() == IloAlgorithm::Infeasible){
                cout << "S.O.S. This should not happen." << endl;
            }
            else if (cplexmTHPMIP.getStatus() == IloAlgorithm::Optimal){
                TPMIP_Obj = cplexmTHPMIP.getObjValue();
                            
                IloNumArray tcysol(env);
                IloNumArray tchsol(env);
                cplexmTHPMIP.getValues(tcysol,tcyvar);
                cplexmTHPMIP.getValues(tchsol,tchvar);

                //Get3rdStageSol(Cycles3rdSol, Chains3rdSol, tcysol, tchsol);
                
                if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
                    if (ThisWork(tcysol, tchsol, vinFirstStage) == true) break;
                }else{//Literature's approach
                    if (Literature(tcysol, tchsol) == true) break;
                }
                
                tcysol.end();
                tchsol.end();
            }else{
                GlobalIte2ndStage += Ite2ndS;
                tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                Print2ndStage("TimeOut");
            }
        }
        else{
            runCGtrue++;
            if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
                if (ThisWork(tcysolColGen, tchsolColGen, vinFirstStage) == true) break;
            }else{//Literature's approach
                if (Literature(tcysolColGen, tchsolColGen) == true) break;
            }
        }

    }
    GlobalIte2ndStage += Ite2ndS;
    tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
    if (LeftTime <= 0) {
        Print2ndStage("TimeOut");
    }
    
    
    //Check whether a new scenario was found
    //if (SPMIP_Obj < FPMIP_Obj){
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
        Print2ndStage("TimeOut");
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
        vector<IndexGrandSubSol>SolFirstStage;
        IloNumArray xsol(env, ListCycles.size());
        cplexRobust.getValues(xsol,X_c);
        for (int i = 0; i < xsol.getSize(); i++){
            if (xsol[i] > 0.9){
                SolFirstStage.push_back(IndexGrandSubSol(ListCycles[i].get_c(), ListCycles[i].get_c().size()));
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
            SolFirstStage.push_back(IndexGrandSubSol(vChains[i], vChains[i].size() - 1));
        }
        if (SPMIP_Obj < FPMIP_Obj){ //Send scenario to 1st. stage
            //Call 2nd. stage
            tStart2ndS = clock();
            Chains2ndStage = Get2ndStageChains (SolFirstStage, RecoursePolicy);
            Cycles2ndStage = Get2ndStageCycles (SolFirstStage, RecoursePolicy);
            tTotalFindingCyCh+= (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            GrandSubProb.end();
            cplexGrandSubP.end();
            Elms2ndPhase.clear();
            Const2ndPhase.clear();
            AllConst2ndPhase.clear();
            Cycles2ndTo3rd.clear();
            Cycles3rdTo2nd.clear();
            Chains2ndTo3rd.clear();
            Chains3rdTo2nd.clear();
            RecoSolCovering.clear();
            RecoTotalWCovering.clear();
            ub_tcyvar.clear();
            ub_tchvar.clear();
            SPMIP_Obj = 0;
            mTHPMIP.end();
            cplexmTHPMIP.end();
            Ite2ndS = 0;
            Ite1stStage++;
            LOWEST_TPMIP_Obj = INT_MAX;
            if (THP_Method != "Benders"){
                GrandSubProbMaster(Cycles2ndStage,Chains2ndStage,SolFirstStage);
            }
            else{
                GrandSubProbMaster2(Cycles2ndStage,Chains2ndStage,SolFirstStage);
            }
        }
        else if (SPMIP_Obj == FPMIP_Obj){
            //Robust Solution found
            Print2ndStage("Optimal");
        }
        else{
            cout << endl << "This should never happen";
        }
    }
    else{
        Print2ndStage("TimeOut");
    }
        
    
    
}
vector<int> SelectConst2ndPh(vector<coverConst>&AllConst2ndPhase, vector<double>&RecoTotalWCovering){
    vector<int>order;
    vector<pair<int,double>>aux;
    for (int i = 0; i < AllConst2ndPhase.size(); i++){
        aux.push_back(make_pair(i, RecoTotalWCovering[i]));
    }
    sort(aux.begin(), aux.end(), sortduals);
    
    for (int i = 0; i < 150; i++){
        order.push_back(aux[i].first);
    }
    return order;
}
bool Problem::ThisWork(IloNumArray& tcysol, IloNumArray& tchsol, vector<int>&vinFirstStage){
        
    double Actual_THPObjective = 0;
    
    for (int i = 0; i < tchsol.getSize(); i++){
        if (tchsol[i] > 0.9){
            //Compute actual THP Objective
            //ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
            map<int,bool>::iterator find0 = ub_tchvar.find(i);
            if (find0 == ub_tchvar.end()){
                Actual_THPObjective+= Chains2ndStage[i].AccumWeight;
            }
        }
    }
    
    for (int i = 0; i < tcysol.getSize(); i++){
        if (tcysol[i] > 0.9){
            //Compute actual THP Objective
            //ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
            map<int,bool>::iterator find0 = ub_tcyvar.find(i);
            if (find0 == ub_tcyvar.end()){
                Actual_THPObjective+= Cycles2ndStage[i].get_Many();
            }
        }
    }
    
    TPMIP_Obj = Actual_THPObjective;
    
    if (SPMIP_Obj == TPMIP_Obj){
        OptFailedArcs = FailedArcs;
        OptFailedVertices = FailedVertices;
        OutforBound++;
        return true; //Upper and Lower bound match
    }
    
    if (TPMIP_Obj < LOWEST_TPMIP_Obj){
        LOWEST_TPMIP_Obj = TPMIP_Obj;
        OptFailedArcs = FailedArcs;
        OptFailedVertices = FailedVertices;
        Get3rdStageSol(Cycles3rdSol, Chains3rdSol, tcysol, tchsol);
        if (Cycles3rdSol.size() == 0 && Chains3rdSol.size() == 0){
            SPMIP_Obj = LOWEST_TPMIP_Obj;
            return true; // No recourse solution found;
        }
        for (int i = 0; i < RecoSolCovering.size(); i++){
            int RHS;
            RHS = Update_RHS_Covering(i);
            if (RHS > AtLeastOneFails[i].getLB()){
                AtLeastOneFails[i].setLB(RHS);
                Const2ndPhase[i].set_RHS(RHS);
            }
        }
    }

    //Add AtLeastOneFails Cut
    if (THP_Method == "Covering"){
        GetAtLeastOneFails(tcysol, tchsol);
    }
    else{
        GetAtLeastOneFailsTwo(tcysol, tchsol);
    }
    
    
    //Resolve 2nd. Phase
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    tStartHeu = clock();
    double ratio = tTotalOptP/Ite1stStage;
    if (Ite2ndS >= 150 && THP_Bound != "NoBound"){
        //Return to optimality
        IteOptP++;
        if (Ite1stStage == 1) IteOptPIte1stis1++;
        Cycles2ndTo3rd.clear();
        Chains2ndTo3rd.clear();
        Cycles3rdTo2nd.clear();
        Chains3rdTo2nd.clear();
        KEPSols2ndStage.clear();
        vector<KEPSol> KEPSol2ndStageUnique;
        KEPSol2ndStageUnique.push_back(KEPSol());
        int counter = -1;
        for (int i = 0; i <  AllConst2ndPhase.size(); i++){
            KEPSols2ndStage.push_back(KEPSol());
            for (int j = 0; j < AllConst2ndPhase[i].get_cycles3rd().size(); j++){
                auto it = Cycles3rdTo2nd.find(AllConst2ndPhase[i].get_cycles3rd()[j]);
                if (it == Cycles3rdTo2nd.end()){
                    counter++;
                    Cycles2ndTo3rd[counter] = AllConst2ndPhase[i].get_cycles3rd()[j];
                    Cycles3rdTo2nd[AllConst2ndPhase[i].get_cycles3rd()[j]] = counter;
                    KEPSol2ndStageUnique.back().cycles.push_back(counter);
                }
                KEPSols2ndStage.back().cycles.push_back(Cycles3rdTo2nd[AllConst2ndPhase[i].get_cycles3rd()[j]]);
            }
        }
        counter = -1;
        for (int i = 0; i <  AllConst2ndPhase.size(); i++){
            for (int j = 0; j < AllConst2ndPhase[i].get_chains3rd().size(); j++){
                auto it = Chains3rdTo2nd.find(AllConst2ndPhase[i].get_chains3rd()[j]);
                if (it == Chains3rdTo2nd.end()){
                    counter++;
                    Chains2ndTo3rd[counter] = AllConst2ndPhase[i].get_chains3rd()[j];
                    Chains3rdTo2nd[AllConst2ndPhase[i].get_chains3rd()[j]] = counter;
                    KEPSol2ndStageUnique.back().chains.push_back(counter);
                }
                KEPSols2ndStage[i].chains.push_back(Chains3rdTo2nd[AllConst2ndPhase[i].get_chains3rd()[j]]);
            }
        }
        //Call GrandSubPromAux
        GrandSubProMastermAux(KEPSols2ndStage, KEPSol2ndStageUnique);
        tTotalOptP+= (clock() - tStartHeu)/double(CLOCKS_PER_SEC);

    }else{
        bool runH = Heuristcs2ndPH();
        tTotalHeu += (clock() - tStartHeu)/double(CLOCKS_PER_SEC);
        if (runH == false){
            tStartMP2ndPH = clock();
            LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
            if (LeftTime < 0){
                GlobalIte2ndStage += Ite2ndS;
                tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                Print2ndStage("TimeOut");
            }
            cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, LeftTime);
            cplexGrandSubP.solve();
            tTotalMP2ndPH += (clock() - tStartMP2ndPH)/double(CLOCKS_PER_SEC);

            if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
                //cout << "2nd. stage solved" << endl;
                SPMIP_Obj = LOWEST_TPMIP_Obj;
                OutforInfeas++;
                return true;
            }
            else if (cplexGrandSubP.getStatus() == IloAlgorithm::Feasible || cplexGrandSubP.getStatus() == IloAlgorithm::Optimal){
                //Get failing arcs and failing vertices
                vertex_sol = IloNumArray(env, Nodes);
                cplexGrandSubP.getValues(vertex_sol,vertex);
                
                arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
                for (int f = 0; f < arc_sol.getSize(); f++){
                    arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
                    cplexGrandSubP.getValues(arc_sol[f],arc[f]);
                }
            }
            else{
                GlobalIte2ndStage += Ite2ndS;
                tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                Print2ndStage("TimeOut");
            }
        }
        else{
            runHeuristicstrue++;
            //Build solution
            vertex_sol = IloNumArray(env, Nodes);
            for (int j = 0; j < scenarioHeuristics.size(); j++){
                if (scenarioHeuristics[j].first == - 1){
                    vertex_sol[scenarioHeuristics[j].second] = 1;
                }
            }
            arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
            for (int f = 0; f < arc_sol.getSize(); f++){
                arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            }
            for (int j = 0; j < scenarioHeuristics.size(); j++){
                if (scenarioHeuristics[j].first != - 1){
                    for (int i = 0; i < AdjacencyList[scenarioHeuristics[j].first].getSize(); i++){
                        if (AdjacencyList[scenarioHeuristics[j].first][i] == scenarioHeuristics[j].second + 1){
                            arc_sol[scenarioHeuristics[j].first][i] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
    
    //Change variables' UB in 3rd phase
    ub_tcyvar.clear(); // <cycle number, 0 or 1>
    ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
    for(int i = 0; i < tcyvar.getSize(); i++){
        map<int,bool>::iterator find = ub_tcyvar.find(i);
        if (find != ub_tcyvar.end()){
            tcyvar[i].setUB(0);
        }
        else{
            tcyvar[i].setUB(1);
        }
    }
    
    ub_tchvar.clear(); // <chain number, 0 or 1>
    ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
    for(int i = 0; i < tchvar.getSize(); i++){
        map<int,bool>::iterator find = ub_tchvar.find(i);
        if (find != ub_tchvar.end()){
            tchvar[i].setUB(0);
        }
        else{
            tchvar[i].setUB(1);
        }
    }
    
    return false;
}
bool Problem::Literature(IloNumArray& tcysol, IloNumArray& tchsol){
    //Create new solution
    
    KEPSols2ndStage.push_back(KEPSol());
   
    vector<KEPSol> KEPSol2ndStageMissing;
    KEPSol2ndStageMissing.push_back(KEPSol());
    
    ub_tcyvar.clear(); // <cycle number, 0 or 1>
    ub_tchvar.clear(); // <cycle number, 0 or 1>
    double Actual_THPObjective = 0;
    
    //Check whether the selected cycles and chains already exist in the 2nd. stage
    for (int i = 0; i < tcysol.getSize(); i++){
        if (tcysol[i] > 0.9){
            //Compute actual THP Objective
            ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
            map<int,bool>::iterator find0 = ub_tcyvar.find(i);
            if (find0 == ub_tcyvar.end()){
                Actual_THPObjective+= Cycles2ndStage[i].get_Many();
            }
            //Check in Cycles3rdTo2nd
            map<int,int>::iterator find = Cycles3rdTo2nd.find(i);
            if (find == Cycles3rdTo2nd.end()){
                //Create new cycle column
                IloNumColumn col(env);
                //Create new cycle variable
                cyvar.add(IloIntVar(col, 0, 1));
                string name = "x." + to_string(int(Cycles3rdTo2nd.size() + 1));
                const char* varName = name.c_str();
                cyvar[cyvar.getSize() - 1].setName(varName);
                cyvar[cyvar.getSize() - 1].setBounds(0, 1);
                //Increase counter
                int aux = int(Cycles3rdTo2nd.size());
                Cycles3rdTo2nd[i] = aux;
                Cycles2ndTo3rd[Cycles3rdTo2nd[i]] = i;
                //Add cycle to KEPSols2ndStage
                KEPSols2ndStage.back().cycles.push_back(Cycles3rdTo2nd[i]);
                KEPSol2ndStageMissing.back().cycles.push_back(Cycles3rdTo2nd[i]);
            }
            else{
                //Add cycle to KEPSols2ndStage
                KEPSols2ndStage.back().cycles.push_back(Cycles3rdTo2nd[i]);
            }
        }
    }

    for (int i = 0; i < tchsol.getSize(); i++){
        if (tchsol[i] > 0.9){
            //Compute actual THP Objective
            ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
            map<int,bool>::iterator find0 = ub_tchvar.find(i);
            if (find0 == ub_tchvar.end()){
                Actual_THPObjective+= Chains2ndStage[i].AccumWeight;
            }
            //Check in Cycles3rdTo2nd
            map<int,int>::iterator find = Chains3rdTo2nd.find(i);
            if (find == Chains3rdTo2nd.end()){
                //Create new chain column
                IloNumColumn col(env);
                //Create new chain variable
                chvar.add(IloIntVar(col, 0, 1));
                string name = "x." + to_string(int(Chains3rdTo2nd.size() + 1));
                const char* varName = name.c_str();
                chvar[chvar.getSize() - 1].setName(varName);
                chvar[chvar.getSize() - 1].setBounds(0, 1);
                //Increase counter
                int aux = int(Chains3rdTo2nd.size());
                Chains3rdTo2nd[i] = aux;
                Chains2ndTo3rd[Chains3rdTo2nd[i]] = i;
                //Add cycle to KEPSols2ndStage
                KEPSols2ndStage.back().chains.push_back(Chains3rdTo2nd[i]);
                KEPSol2ndStageMissing.back().chains.push_back(Chains3rdTo2nd[i]);
            }
            else{
                //Add cycle to KEPSols2ndStage
                KEPSols2ndStage.back().chains.push_back(Chains3rdTo2nd[i]);
            }
        }
    }
    
    
    //Obtain true TPMIP_Obj if THP_Method == Literature
    TPMIP_Obj = Actual_THPObjective;
    
    //Create Bounding Constraint
    Const11b(KEPSols2ndStage);
    
    GrandSubProb.add(vBoundConstraint[vBoundConstraint.getSize() - 1]);
    exprBound.end();
    
    
    //Create Active Cols Constraint
    Const11c(KEPSol2ndStageMissing);
    GrandSubProb.add(vFailedMatches[vFailedMatches.getSize() - 1]);
    exprVxtArcsCH.end();
    exprVxtArcsCY.end();
    
    //Resolve 2ndPhase
    //cplexGrandSubP.exportModel("GrandSubP2.lp");
    cplexGrandSubP.solve();
    if (cplexGrandSubP.isPrimalFeasible() == false){
        cout << "IT SHOULD NEVER HAPPEN" << endl;
    }
    else{
        SPMIP_Obj = cplexGrandSubP.getObjValue();
        SPMIP_Obj = round(SPMIP_Obj);
        if (SPMIP_Obj >= TPMIP_Obj){
            OptFailedArcs = FailedArcs;
            OptFailedVertices = FailedVertices;
            return true;
        }

        IloNumArray cyvar_sol(env, Cycles2ndTo3rd.size());
        IloNumArray chvar_sol(env, Chains2ndTo3rd.size());
        cplexGrandSubP.getValues(cyvar_sol,cyvar);
        cplexGrandSubP.getValues(chvar_sol,chvar);
        
//        cout << "Sol. 2nd. Phase: " << endl;
//        for (int i = 0; i < cyvar_sol.getSize(); i++){
//            if (cyvar_sol[i] > 0.1) cout << cyvar[i].getName() << endl;
//        }
//        for (int i = 0; i < chvar_sol.getSize(); i++){
//            if (chvar_sol[i] > 0.1) cout << chvar[i].getName() << endl;
//        }
        cyvar_sol.end();
        chvar_sol.end();
        //Get failing arcs and failing vertices
        vertex_sol = IloNumArray(env, Nodes);
        cplexGrandSubP.getValues(vertex_sol,vertex);
        
        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            cplexGrandSubP.getValues(arc_sol[f],arc[f]);
        }
        GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
        
        //Change variables' UB in 3rd phase
        
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
        
        ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
        for(int i = 0; i < tchvar.getSize(); i++){
            map<int,bool>::iterator find = ub_tchvar.find(i);
            if (find != ub_tchvar.end()){
                //Change ObjCoefficients
                ObjTHP.setLinearCoef(tchvar[i], 1);
            }
            else{
                //Change ObjCoefficients
                int n = int(Chains2ndStage[i].AccumWeight);
                ObjTHP.setLinearCoef(tchvar[i], (n*Nodes + 1));
            }
        }
        mTHPMIP.add(ObjTHP);
    
    }
    return false;
    
}
void Problem::GrandSubProMastermAux(vector<KEPSol>KEPSols, vector<KEPSol>KEPUniqueEx){
    // Create model
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    if (LeftTime < 0){
        GlobalIte2ndStage += Ite2ndS;
        tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
        Print2ndStage("TimeOut");
    }
    IloModel GrandSubProbAux(env);
    IloCplex cplexGrandAux(GrandSubProbAux);
    cplexGrandAux.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexGrandAux.setParam(IloCplex::Param::Threads, 1);
    cplexGrandAux.setOut(env.getNullStream());
    
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
    NumVar2D Auxarc(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        Auxarc[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, 1, ILOINT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            mapArcs[make_pair(i, AdjacencyList[i][j] - 1)] = j;
            SetName2(arc[i][j], "arc", i + 1, AdjacencyList[i][j]);
        }
    }
    
    //Create vertex variables
    IloNumVarArray Auxvertex(env, Nodes, 0, 1, ILOINT);
    for (int i = 0; i < Nodes; i++){
        SetName1Index( Auxvertex[i], "vertex",i + 1);
    }
    
    //Create bounding variable
    Beta = IloNumVar(env);
    SetName1Index(Beta, "Beta", 0);
    
    //Create Bounding Constraint
    vBoundConstraint = IloRangeArray(env);
    for (int s = 0; s < KEPSols.size(); s++){
        //Const11b
        exprBound = IloExpr(env,0);
        for (int i = 0; i < KEPSols[s].cycles.size(); i++){
            int w = Cycles2ndStage[Cycles2ndTo3rd[KEPSols[s].cycles[i]]].get_Many();
            exprBound+= w*cyvar[KEPSols[s].cycles[i]];
        }
        for (int i = 0; i < KEPSols[s].chains.size(); i++){
            int w = Chains2ndStage[Chains2ndTo3rd[KEPSols[s].chains[i]]].AccumWeight;
            exprBound+= w*chvar[KEPSols[s].chains[i]];
        }
        //cout << Beta.getName() << ">=" << exprBound << endl;
        string name = "Const11b."  + to_string(vBoundConstraint.getSize());
        const char* cName = name.c_str();
        vBoundConstraint.add(IloRange(env, 0, Beta - exprBound , IloInfinity, cName));
    }
    GrandSubProbAux.add(vBoundConstraint);
    exprBound.end();
    
    //Create Active Cols Constraint
    vFailedMatches = IloRangeArray(env);
    for (int i = 0; i < KEPUniqueEx.back().cycles.size(); i++){
        exprVxtArcsCY = IloExpr (env, 0);
        int u = KEPUniqueEx.back().cycles[i];
        for (int j = 0; j < Cycles2ndStage[Cycles2ndTo3rd[u]].get_c().size(); j++){
            exprVxtArcsCY += Auxvertex[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]];
            if (j <= Cycles2ndStage[Cycles2ndTo3rd[u]].get_c().size() - 2){
                int v = mapArcs[make_pair(Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j], Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j + 1])];
                exprVxtArcsCY += Auxarc[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]][v];
            }
            else{
                int v = mapArcs[make_pair(Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j], Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[0])];
                exprVxtArcsCY += Auxarc[Cycles2ndStage[Cycles2ndTo3rd[u]].get_c()[j]][v];
            }
        }
        //cout << cyvar[u] + exprVxtArcsCY << ">=" << 1 << endl;
        string name = "Const11c.CY"  + to_string(vFailedMatches.getSize());
        const char* cName = name.c_str();
        vFailedMatches.add(IloRange(env, 1, cyvar[u] + exprVxtArcsCY , IloInfinity, cName));
    }
    for (int i = 0; i < KEPUniqueEx.back().chains.size(); i++){
        exprVxtArcsCH = IloExpr (env, 0);
        int u = KEPUniqueEx.back().chains[i];
        for (int j = 0; j < Chains2ndStage[Chains2ndTo3rd[u]].Vnodes.size(); j++){
            exprVxtArcsCH += Auxvertex[Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex];
            if (j <= Chains2ndStage[Chains2ndTo3rd[u]].Vnodes.size() - 2){
                int v = mapArcs[make_pair(Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex, Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j + 1].vertex)];
                exprVxtArcsCH += Auxarc[Chains2ndStage[Chains2ndTo3rd[u]].Vnodes[j].vertex][v];
            }
        }
        //cout << chvar[u] + exprVxtArcsCH << ">=" << 1 << endl;
        string name = "Const11c.CH"  + to_string(vFailedMatches.getSize());
        const char* cName = name.c_str();
        vFailedMatches.add(IloRange(env, 1, chvar[u] + exprVxtArcsCH , IloInfinity, cName));
//        GrandSubProb.add(vFailedMatches[vFailedMatches.getSize() - 1]);
//        exprVxtArcsCH.end();
    }
    GrandSubProbAux.add(vFailedMatches);
    exprVxtArcsCH.end();
    exprVxtArcsCY.end();
    
    //Create arcs and vertex constraints
    IloExpr sumVertices (env, 0);
    for (int i = 0; i < Nodes; i++){
        sumVertices+=  Auxvertex[i];
    }
    string name;
    const char* cName;
    name = "VtxSum";
    cName = name.c_str();
    GrandSubProbAux.add(IloRange(env, -IloInfinity, sumVertices, MaxVertexFailures, cName));
    
    IloExpr sumArcs (env, 0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumArcs+= Auxarc[i][j];
        }
    }
    name = "ArcSum";
    cName = name.c_str();
    GrandSubProbAux.add(IloRange(env, -IloInfinity, sumArcs, MaxArcFailures, cName));
    
    //Add objective
    name = "Obj_GrandSubP";
    GrandSubProbAux.add(IloObjective(env, Beta, IloObjective::Minimize, name.c_str()));
    
    //cplexGrandSubP.exportModel("GrandSubP2.lp");
    cplexGrandAux.solve();
    if (cplexGrandAux.getStatus() == IloAlgorithm::Infeasible){
        cout << "S.O.S. This should not happen." << endl;
    }
    else if (cplexGrandAux.getStatus() == IloAlgorithm::Optimal){
        
        SPMIP_Obj = cplexGrandAux.getValue(Beta);
        
        //Retrieve solution
        IloNumArray cyvar_sol(env, Cycles2ndTo3rd.size());
        IloNumArray chvar_sol(env, Chains2ndTo3rd.size());
        cplexGrandAux.getValues(cyvar_sol,cyvar);
        cplexGrandAux.getValues(chvar_sol,chvar);
        
        vertex_sol = IloNumArray(env, Nodes);
        cplexGrandAux.getValues(vertex_sol,Auxvertex);
        
        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            cplexGrandAux.getValues(arc_sol[f],Auxarc[f]);
        }
    }
    else{
        GlobalIte2ndStage += Ite2ndS;
        tTotal2ndS += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
        Print2ndStage("TimeOut");
    }
    GrandSubProbAux.end();
    cplexGrandAux.end();
}
void Problem::GetNewBetaCut(IloNum TPMIP_Obj, map<pair<int,int>, bool> FailedArcs, map<int, bool> FailedVertices){
    IloExpr expr(env, 0);
    for (auto it = mapArcs.begin(); it != mapArcs.end(); it++){
        for (auto it2 = FailedArcs.begin(); it2 != FailedArcs.end(); it2++){
            auto itf = FailedArcs.find(it->first);
            if (itf == FailedArcs.end()){
                expr+= arc[it->first.first][it->second];
            }
        }
    }
    for (int i = 0; i < Nodes; i++){
        for (auto it2 = FailedVertices.begin(); it2 != FailedVertices.end(); it2++){
            auto itf = FailedVertices.find(i);
            if (itf == FailedVertices.end()){
                expr+= vertex[i];
            }
        }
    }
    string name = "Beta." + to_string(ConsBeta.getSize() + 1);
    cout << Beta + expr << endl;
    ConsBeta.add(IloRange(env, TPMIP_Obj, Beta + TPMIP_Obj*expr, IloInfinity, name.c_str()));
    GrandSubProb.add(ConsBeta[ConsBeta.getSize() - 1]);
//    for (int i = 0; i < tcysol.getSize(); i++){
//        if (tcysol[i] > 0.9){
//            expr-= Cycles2ndStage[i].get_Many()*cyvar[Cycles3rdTo2nd[i]];
//        }
//    }
//    for (int i = 0; i < tchsol.getSize(); i++){
//        if (tchsol[i] > 0.9){
//            expr-= Chains2ndStage[i].AccumWeight*chvar[Chains3rdTo2nd[i]];
//        }
//    }
//    string name = "Beta." + to_string(ConsBeta.getSize() + 1);
//    cout << Beta + expr << endl;
//    ConsBeta.add(IloRange(env, 0, Beta + expr, IloInfinity, name.c_str()));
//    GrandSubProb.add(ConsBeta[ConsBeta.getSize() - 1]);
}
void Problem::GetAtLeastOneFails(IloNumArray& tcysol, IloNumArray& tchsol){

    IloExpr expr(env);
    Cycles3rdSol.clear();
    Chains3rdSol.clear();
    vector<int>cycles;
    vector<int>chains;
    vector<int>Allcycles;
    vector<int>Allchains;
    for (int i = 0; i < tcysol.getSize(); i++){
        map<int,bool>::iterator find0 = ub_tcyvar.find(i);
        if (tcysol[i] > 0.9){
            Allcycles.push_back(i);
        }
        if (tcysol[i] > 0.9 && find0 == ub_tcyvar.end()){
            cycles.push_back(i);
            Cycles3rdSol.push_back(Cycles2ndStage[i]);
            for (int j = 0; j < Cycles2ndStage[i].get_c().size(); j++){
                int u = Cycles2ndStage[i].get_c()[j];
                auto it = Elms2ndPhase.find(make_pair(-1,u));
                if (it == Elms2ndPhase.end()){
                    Elms2ndPhase[make_pair(-1,u)] = coveringElements();
                }
                Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
                Elms2ndPhase[make_pair(-1,u)].add_weight(Cycles2ndStage[i].get_Many());
                Elms2ndPhase[make_pair(-1,u)].add_map(int(AtLeastOneFails.getSize()),i);
                expr+= vertex[u];
                if (j == Cycles2ndStage[i].get_c().size() - 1){
                    int s = mapArcs[make_pair(u,Cycles2ndStage[i].get_c()[0])];
                    expr+= arc[u][s];
                    auto it = Elms2ndPhase.find(make_pair(u,Cycles2ndStage[i].get_c()[0]));
                    if (it == Elms2ndPhase.end()){
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])] = coveringElements();
                    }
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_const(int(AtLeastOneFails.getSize()));
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_weight(Cycles2ndStage[i].get_Many());
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_map(int(AtLeastOneFails.getSize()), i);
                }
                else{
                    int s = mapArcs[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])];
                    expr+= arc[u][s];
                    auto it = Elms2ndPhase.find(make_pair(u,Cycles2ndStage[i].get_c()[j + 1]));
                    if (it == Elms2ndPhase.end()){
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])] = coveringElements();
                    }
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_const(int(AtLeastOneFails.getSize()));
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_weight(Cycles2ndStage[i].get_Many());
                    Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_map(int(AtLeastOneFails.getSize()), i);
                }
            }
        }
    }
    for (int i = 0; i < tchsol.getSize(); i++){
        map<int,bool>::iterator find0 = ub_tchvar.find(i);
        if (tchsol[i] > 0.9){
            Allchains.push_back(i);
        }
        if (tchsol[i] > 0.9 && find0 == ub_tchvar.end()){
            chains.push_back(i);
            Chains3rdSol.push_back(Chains2ndStage[i]);
            for (int j = 0; j < Chains2ndStage[i].Vnodes.size(); j++){
                int u = Chains2ndStage[i].Vnodes[j].vertex;
                auto it = Elms2ndPhase.find(make_pair(-1,u));
                if (it == Elms2ndPhase.end()){
                    Elms2ndPhase[make_pair(-1,u)] = coveringElements();
                }
                Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
                Elms2ndPhase[make_pair(-1,u)].add_weight(Chains2ndStage[i].AccumWeight);
                Elms2ndPhase[make_pair(-1,u)].add_map(int(AtLeastOneFails.getSize()), i);
                expr+= vertex[u];
                if (j <= Chains2ndStage[i].Vnodes.size() - 2){
                    int s = mapArcs[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)];
                    expr+= arc[u][s];
                    auto it = Elms2ndPhase.find(make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex));
                    if (it == Elms2ndPhase.end()){
                        Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)] = coveringElements();
                    }
                    Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_const(int(AtLeastOneFails.getSize()));
                    Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_weight(Chains2ndStage[i].AccumWeight);
                    Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_map(int(AtLeastOneFails.getSize()), i);
                }
            }
        }
    }
    if (Cycles3rdSol.size() > 0 || Chains3rdSol.size() > 0){
        int RHS = 1;
        if (Ite2ndS >= 1){
            //Sort Cycles3rdSol and Chains3rdSol
            RecoSolCovering.push_back(vector<double>());
            double TotalW = 0;
            for (int i = 0; i < Cycles3rdSol.size(); i++){
                RecoSolCovering.back().push_back(Cycles3rdSol[i].get_Many());
                TotalW+= Cycles3rdSol[i].get_Many();
            }
            for (int i = 0; i < Chains3rdSol.size(); i++){
                RecoSolCovering.back().push_back(Chains3rdSol[i].AccumWeight);
                TotalW+= Chains3rdSol[i].AccumWeight;
            }
            RecoTotalWCovering.push_back(TotalW);
            sort(RecoSolCovering.back().begin(), RecoSolCovering.back().end(), sortdouble);
            
        }

        if (RecoTotalWCovering.back() >= LOWEST_TPMIP_Obj){
            RHS = Update_RHS_Covering(int (RecoSolCovering.size() - 1));
        }
        //Add new constraint to Const2ndPhase
        Const2ndPhase.push_back(coverConst(cycles, chains));
        Const2ndPhase.back().set_RHS(RHS);
        AllConst2ndPhase.push_back(coverConst(Allcycles, Allchains));
        
        string name = "AtLeastOneFails_" + to_string(AtLeastOneFails.getSize() + 1);
        const char* cName = name.c_str();
        //cout << expr << endl;
        AtLeastOneFails.add(IloRange(env, RHS, expr, IloInfinity, cName));
        GrandSubProb.add(AtLeastOneFails[AtLeastOneFails.getSize() - 1]);
        expr.end();
    }
    
}
void Problem::GetAtLeastOneFailsTwo(IloNumArray& tcysol, IloNumArray& tchsol){

int counter = 0;
while(counter <= 2){
    counter++;
    IloExpr expr(env);
    Cycles3rdSol.clear();
    Chains3rdSol.clear();
    vector<int>cycles;
    vector<int>chains;
    vector<int>Allcycles;
    vector<int>Allchains;
    for (int i = 0; i < tcysol.getSize(); i++){
        map<int,bool>::iterator find0 = ub_tcyvar.find(i);
        if (tcysol[i] > 0.9 && counter == 1){
            Allcycles.push_back(i);
        }
        if (tcysol[i] > 0.9){
            bool skip = false;
            if (counter >= 2 && find0 != ub_tcyvar.end()) skip = true;
            if (skip == false){
                cycles.push_back(i);
                Cycles3rdSol.push_back(Cycles2ndStage[i]);
                for (int j = 0; j < Cycles2ndStage[i].get_c().size(); j++){
                    int u = Cycles2ndStage[i].get_c()[j];
                    auto it = Elms2ndPhase.find(make_pair(-1,u));
                    if (it == Elms2ndPhase.end()){
                        Elms2ndPhase[make_pair(-1,u)] = coveringElements();
                    }
                    Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
                    Elms2ndPhase[make_pair(-1,u)].add_weight(Cycles2ndStage[i].get_Many());
                    Elms2ndPhase[make_pair(-1,u)].add_map(int(AtLeastOneFails.getSize()),i);
                    expr+= vertex[u];
                    if (j == Cycles2ndStage[i].get_c().size() - 1){
                        int s = mapArcs[make_pair(u,Cycles2ndStage[i].get_c()[0])];
                        expr+= arc[u][s];
                        auto it = Elms2ndPhase.find(make_pair(u,Cycles2ndStage[i].get_c()[0]));
                        if (it == Elms2ndPhase.end()){
                            Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])] = coveringElements();
                        }
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_const(int(AtLeastOneFails.getSize()));
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_weight(Cycles2ndStage[i].get_Many());
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[0])].add_map(int(AtLeastOneFails.getSize()), i);
                    }
                    else{
                        int s = mapArcs[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])];
                        expr+= arc[u][s];
                        auto it = Elms2ndPhase.find(make_pair(u,Cycles2ndStage[i].get_c()[j + 1]));
                        if (it == Elms2ndPhase.end()){
                            Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])] = coveringElements();
                        }
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_const(int(AtLeastOneFails.getSize()));
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_weight(Cycles2ndStage[i].get_Many());
                        Elms2ndPhase[make_pair(u,Cycles2ndStage[i].get_c()[j + 1])].add_map(int(AtLeastOneFails.getSize()), i);
                    }
                }
            }
        }
    }
    for (int i = 0; i < tchsol.getSize(); i++){
        map<int,bool>::iterator find0 = ub_tchvar.find(i);
        if (tchsol[i] > 0.9 && counter == 1){
            Allchains.push_back(i);
        }
        if (tchsol[i] > 0.9){
            bool skip = false;
            if (counter >= 2 && find0 != ub_tchvar.end()) skip = true;
            if (skip == false){
                chains.push_back(i);
                Chains3rdSol.push_back(Chains2ndStage[i]);
                for (int j = 0; j < Chains2ndStage[i].Vnodes.size(); j++){
                    int u = Chains2ndStage[i].Vnodes[j].vertex;
                    auto it = Elms2ndPhase.find(make_pair(-1,u));
                    if (it == Elms2ndPhase.end()){
                        Elms2ndPhase[make_pair(-1,u)] = coveringElements();
                    }
                    Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
                    Elms2ndPhase[make_pair(-1,u)].add_weight(Chains2ndStage[i].AccumWeight);
                    Elms2ndPhase[make_pair(-1,u)].add_map(int(AtLeastOneFails.getSize()), i);
                    expr+= vertex[u];
                    if (j <= Chains2ndStage[i].Vnodes.size() - 2){
                        int s = mapArcs[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)];
                        expr+= arc[u][s];
                        auto it = Elms2ndPhase.find(make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex));
                        if (it == Elms2ndPhase.end()){
                            Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)] = coveringElements();
                        }
                        Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_const(int(AtLeastOneFails.getSize()));
                        Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_weight(Chains2ndStage[i].AccumWeight);
                        Elms2ndPhase[make_pair(u,Chains2ndStage[i].Vnodes[j + 1].vertex)].add_map(int(AtLeastOneFails.getSize()), i);
                    }
                }
            }
        }
    }
    if (Cycles3rdSol.size() > 0 || Chains3rdSol.size() > 0){
        int RHS = 1;
        if (Ite2ndS >= 1){
            //Sort Cycles3rdSol and Chains3rdSol
            RecoSolCovering.push_back(vector<double>());
            double TotalW = 0;
            for (int i = 0; i < Cycles3rdSol.size(); i++){
                RecoSolCovering.back().push_back(Cycles3rdSol[i].get_Many());
                TotalW+= Cycles3rdSol[i].get_Many();
            }
            for (int i = 0; i < Chains3rdSol.size(); i++){
                RecoSolCovering.back().push_back(Chains3rdSol[i].AccumWeight);
                TotalW+= Chains3rdSol[i].AccumWeight;
            }
            RecoTotalWCovering.push_back(TotalW);
            sort(RecoSolCovering.back().begin(), RecoSolCovering.back().end(), sortdouble);
            
        }

        if (RecoTotalWCovering.back() >= LOWEST_TPMIP_Obj){
            RHS = Update_RHS_Covering(int (RecoSolCovering.size() - 1));
        }
        //Add new constraint to Const2ndPhase
        Const2ndPhase.push_back(coverConst(cycles, chains));
        Const2ndPhase.back().set_RHS(RHS);
        if(counter == 1) AllConst2ndPhase.push_back(coverConst(Allcycles, Allchains));
        
        string name = "AtLeastOneFails_" + to_string(AtLeastOneFails.getSize() + 1);
        const char* cName = name.c_str();
        //cout << expr << endl;
        AtLeastOneFails.add(IloRange(env, RHS, expr, IloInfinity, cName));
        GrandSubProb.add(AtLeastOneFails[AtLeastOneFails.getSize() - 1]);
        expr.end();
        counter++;
        if (Ite2ndS == 0) break;
    }
}
        
    
}

bool sortdouble(double& c1, double& c2){
    return (c1 > c2);
}
bool sortduals(pair<int,double>& c1, pair<int,double>& c2){
    return (c1.second < c2.second);
}
bool sortCovElms(pair< pair<int,int>, double>& c1, pair< pair<int,int>, double>& c2){
    return (c1.second > c2.second);
}
bool sortPairs(pair<int,int>& c1, pair<int,int>& c2){
    return (c1.second > c2.second);
}
int Problem::Update_RHS_Covering(int row){
    int i = 0, RHS = 0, accum = 0;

    while (true){
        if (RecoTotalWCovering[row] - RecoSolCovering[row][i] - accum >= LOWEST_TPMIP_Obj){
            accum+= RecoSolCovering[row][i];
            i++;
        }
        else{
            RHS = i + 1;
            break;
        }
    }
    
    return RHS;
}
void Problem::Get3rdStageSol(vector<Cycles>&Cycles3rdSol, vector<Chain>&Chains3rdSol, IloNumArray& cyvar_sol3rd, IloNumArray& chvar_sol3rd){
    Cycles3rdSol.clear();
    Chains3rdSol.clear();
    for (int i = 0; i < cyvar_sol3rd.getSize(); i++){
        if (cyvar_sol3rd[i] > 0.9){
            Cycles3rdSol.push_back(Cycles2ndStage[i]);
        }
    }
    for (int i = 0; i < chvar_sol3rd.getSize(); i++){
        if (chvar_sol3rd[i] > 0.9){
            Chains3rdSol.push_back(Chains2ndStage[i]);
        }
    }
}
IloExpr Problem::GetObjTPH(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, string& TPH_Method){
    IloExpr expr(env,0);
    if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
        for(int i = 0; i < Cycles2ndStage.size(); i++){
            int n = int(Cycles2ndStage[i].get_Many());
                expr += n*tcyvar[i];
        }
        for(int i = 0; i < Chains2ndStage.size(); i++){
            int n = int(Chains2ndStage[i].AccumWeight);
                expr += n*tchvar[i];
        }
        return expr;
    }else if (THP_Method == "Benders"){
        map<int,bool>ub_tcyvar; // <cycle number, 0 or 1>
        ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
        for(int i = 0; i < Cycles2ndStage.size(); i++){
            map<int,bool>::iterator find = ub_tcyvar.find(i);
            if (find != ub_tcyvar.end()){
                expr += tcyvar[i];
            }
            else{
                int n = int(Cycles2ndStage[i].get_Many());
                    expr += (n*Nodes + 1)*tcyvar[i];
            }
        }
        map<int,bool>ub_tchvar; // <cycle number, 0 or 1>
        ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
        for(int i = 0; i < Chains2ndStage.size(); i++){
            map<int,bool>::iterator find = ub_tchvar.find(i);
            if (find != ub_tchvar.end()){
                expr += tchvar[i];
            }
            else{
                int n = int(Chains2ndStage[i].AccumWeight);
                expr += (n*Nodes + 1)*tchvar[i];
            }
        }
        //cout << endl << expr << endl;
        return expr;
    }
    else{
        return expr;
    }

}
void Problem::GetScenario(IloNumArray2& arc_sol, IloNumArray& vertex_sol){
    FailedArcs.clear();
    FailedVertices.clear();
    
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        if (vertex_sol[i] > 0.9) FailedVertices[i] = true;
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            if (arc_sol[i][j] > 0.9)  FailedArcs[make_pair(i, AdjacencyList[i][j] - 1)] = true;
        }
    }
    
}
IloNumVarArray Problem::Create_tcyvar(const char* prefix, vector<Cycles>&Cycles2ndStage){
    //Get upper bound of cycles according to scenario
    map<int,bool>ub_tcyvar; // <cycle number, 0 or 1>
    if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
        ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
    }
    
    IloNumVarArray var(env, Cycles2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
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
        else{
            var[i] = IloNumVar(env, 0, 1, ILOINT);
            SetName(var[i], prefix, i + 1);
        }
    }
    return var;
}
map<int,bool> Problem::GetUB_tcyvar(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices){//returns a map of cycles/chains that fail under a scenario
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
    if (THP_Method != "Benders"){
        ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
    }
    
    IloNumVarArray var(env, Chains2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        if (THP_Method != "Benders"){
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
        //cout << expr << endl;
        DisjointTHPArray.add(IloRange(env, -IloInfinity, expr, 1));
    }
    return DisjointTHPArray;
}
void Problem::AddNewCols3rdTo2nd (IloNumArray tcysol, IloNumArray tchsol, map<int,int>& Cycles2ndTo3rd, map<int,int>& Chains2ndTo3rd, map<int,int>& Cycles3rdTo2nd, map<int,int>& Chains3rdTo2nd, int& Cyclenewrow2ndPH, int& Chainnewrow2ndPH, vector<int>& tcysol3rd, vector<int>& tchsol3rd, vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage){
    Cyclenewrow2ndPH = int(Cycles2ndTo3rd.size());
    Chainnewrow2ndPH = int(Chains2ndTo3rd.size());
    tcysol3rd.clear();
    tchsol3rd.clear();
    for (int i = 0; i < tcysol.getSize(); i++){
        if (tcysol[i] > 0.9){
            map<int,int>::iterator it = Cycles3rdTo2nd.find(i);
            if (it == Cycles3rdTo2nd.end()){
                tcysol3rd.push_back(i);
                Cycles2ndTo3rd[int(Cycles2ndTo3rd.size())] = i;
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
            map<int,int>::iterator it = Chains3rdTo2nd.find(i);
            if (it == Chains3rdTo2nd.end()){
                tchsol3rd.push_back(i);
                Chains2ndTo3rd[int(Chains2ndTo3rd.size())] = i;
                Chains3rdTo2nd[i] = int(Chains2ndTo3rd.size() - 1);
                //Update ChainNodeSPH
                for (int j = 0; j < Chains2ndStage[i].Vnodes.size(); j++){
                    ChainNodeSPH[Chains2ndStage[i].Vnodes[j].vertex].push_back(int(Chains2ndTo3rd.size() - 1));
                }
            }
        }
    }
    
}
//////////////////////////// SVIs ////////////////////////////
bool NeededElement(vector<int>&uncoveredc, coveringElements& Ele){
    if (Ele.get_state() == false){
        for (int i = 0; i < Ele.get_coveredconsts().size(); i++){
            for (int j = 0; j < uncoveredc.size(); j++){
                if (Ele.get_coveredconsts()[i] == uncoveredc[j]){
                    return true;
                }
            }
        }
    }
    
    return false;
}
bool Problem::Heuristcs2ndPH(){
    bool another = false, keepgoing = true, complete = false, non2gether = false;
    int UsedArcs = 0;
    int UsedVertices = 0, ele, ite = 0;
    vector<pair<int,int>>tabuList;
    while(ite <= 1 && keepgoing == true){
        scenarioHeuristics.clear();
        tabuList.clear();
        for (auto it = Elms2ndPhase.begin(); it != Elms2ndPhase.end(); it++) it->second.set_state(false);
        for (int i = 0; i < Const2ndPhase.size(); i++){
            Const2ndPhase[i].deleteEls();
        }
        map<pair<int,int>, coveringElements>Eleaux = Elms2ndPhase;
        UsedArcs = 0;
        UsedVertices = 0;
        ite++;
        while (true){
            vector<pair< pair<int,int>, double>>auxCov;
            for (auto it = Eleaux.begin(); it != Eleaux.end(); it++){
                double w = 0;
                if (it->first.first == -1 && UsedVertices < MaxVertexFailures && Eleaux[it->first].get_state() == false){
                    //w = 0.5*(CycleNodeTPH[it->first.second].size() + ChainNodeTPH[it->first.second].size()) +  0.5*it->second.get_coversize();
                    //w = it->second.get_coversize();0.05*(PredMap[it->first.second].size() + AdjacencyList[it->first.second].getSize()) +
                    //w = it->second.get_coversize(); //Total degree
                    w = it->second.get_coversize();
                    auxCov.push_back(make_pair(it->first, w));
                }
                if (it->first.first != -1 && UsedArcs < MaxArcFailures && Eleaux[it->first].get_state() == false){
                    //w = it->second.get_maxw();
                    w = 0.4*(ArcsinChCyTHP[it->first].size()/(Cycles2ndStage.size() + Chains2ndStage.size())) +  0.6*it->second.get_coversize(); //Number Arcs in cycles/chains 2ndStage
                    auxCov.push_back(make_pair(it->first, w));
                }
            }
            //If no Elements to pick from
            if (auxCov.size() == 0){
                break;
            }
            //Sort elements
            sort(auxCov.begin(), auxCov.end(), sortCovElms);
            while(true){
                double search = auxCov[0].second;
                int UpTo = auxCov.size();//In case all values are repeated equally
                another = false;
                //Choose among the most repeated ones
                for (int i = 0; i < auxCov.size(); i++){
                    if (auxCov[i].second != search){
                        UpTo = i;
                        break;
                    }
                }
                ele = rand()%UpTo;
                //Verify ele is not in tabuList
                 if (auxCov[ele].first.first != -1){// It is an arc
                    for (int i = 0; i < tabuList.size(); i++){
                        if (auxCov[ele].first == tabuList[i]){
                            another = true;
                            //remove element from auxCov
                            auxCov.erase(auxCov.begin() + ele);
                            break;
                        }
                    }
                }
                if (another == false){
                    //Update tabuList
                    if (auxCov[ele].first.first == -1){
                        for (int i = 0; i < PredMap[auxCov[ele].first.second].size(); i++){
                            int val = PredMap[auxCov[ele].first.second][i].first;
                            tabuList.push_back(make_pair(val, auxCov[ele].first.second));
                        }
                        for (int i = 0; i < AdjacencyList[auxCov[ele].first.second].getSize(); i++){
                            tabuList.push_back(make_pair(auxCov[ele].first.second, AdjacencyList[auxCov[ele].first.second][i] - 1));
                        }
                    }
                    break;
                }
            }
            //Check whether we should add element to scenario
            bool ans = true;
            if (UsedArcs + UsedVertices >= 1 && auxCov[ele].first.second < Pairs && (UsedArcs >= 1 || UsedVertices >= 2)){
                vector<int> failedVxt, vinFirstStage;
                vector<pair<int,int>>failedArcs;
                for (int i = 0; i < scenarioHeuristics.size();i++){
                    if(scenarioHeuristics[i].first == -1){
                        failedVxt.push_back(scenarioHeuristics[i].second);
                    }
                    else{
                        failedArcs.push_back(scenarioHeuristics[i]);
                    }
                }
                ans = UnMVtxdueToVtx(failedVxt, failedArcs, vinFirstStage, auxCov[ele].first);
            }
            if (ans == true){
                // Add element to scenario
                scenarioHeuristics.push_back(auxCov[ele].first);
                //Seize the element
                Eleaux[auxCov[ele].first].set_state(true);
                if (auxCov[ele].first.first == -1){//It is a vertex
                    UsedVertices++;
                }
                else{
                    UsedArcs++;
                }
                //Cover constraints
                for (int i = 0; i < Eleaux[auxCov[ele].first].get_coveredconsts().size(); i++){
                    Const2ndPhase[Eleaux[auxCov[ele].first].get_coveredconsts()[i]].add_cover(auxCov[ele].first);
                }
            }
            else{
                Eleaux[auxCov[ele].first].set_state(true);//true so that it is not selected
                non2gether = true;
            }
            keepgoing = false;
            //Check whether there is an unsatisfied constraint
            for (int i = 0; i < Const2ndPhase.size(); i++){
                if (Const2ndPhase[i].get_coversize() < Const2ndPhase[i].get_RHS()){
                    keepgoing = true;
                }
            }
            if (keepgoing == false) {
                complete = true;
                //Check weather we have used up the failure budget. If not, complete
                if (UsedVertices == MaxVertexFailures && UsedArcs == MaxArcFailures){
                    break;
                }
                else{
                    if (auxCov.size() == 0) break;
                    keepgoing = true;
                }
            }
        }
    }
    
    return complete;
}
void Problem::HeuristicsStart2ndPH(map<int,int>& Cycles2ndTo3rd, map<int,int>& Chains2ndTo3rd, vector<int>&ListSelVex){
    int UsedArcs = 0;
    int UsedVertices = 0;
    vector<pair< pair<int,int>, double>>AllArcs;
    vector<pair<int,int>>Arcsin1stStageSol;
    vector<pair<int, int>>AllVex;
    
    if (MaxArcFailures > 0){
        for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
            for (int j = 0; j < Cycles2ndStage[Cycles2ndTo3rd[i]].get_c().size(); j++){
                int u = Cycles2ndStage[Cycles2ndTo3rd[i]].get_c()[j];
                if (j < Cycles2ndStage[Cycles2ndTo3rd[i]].get_c().size() - 1){
                    int v = Cycles2ndStage[Cycles2ndTo3rd[i]].get_c()[j + 1];
                    Arcsin1stStageSol.push_back(make_pair(u,v));
                }
                else{
                    int v = Cycles2ndStage[Cycles2ndTo3rd[i]].get_c()[0];
                    Arcsin1stStageSol.push_back(make_pair(u,v));
                }
            }
        }
        for (int i = 0; i < Chains2ndTo3rd.size(); i++){
            for (int j = 0; j < Chains2ndStage[Chains2ndTo3rd[i]].Vnodes.size(); j++){
                int u = Chains2ndStage[Chains2ndTo3rd[i]].Vnodes[j].vertex;
                if (j < Chains2ndStage[Chains2ndTo3rd[i]].Vnodes.size() - 1){
                    int v = Chains2ndStage[Chains2ndTo3rd[i]].Vnodes[j + 1].vertex;
                    Arcsin1stStageSol.push_back(make_pair(u,v));
                }
            }
        }
        for (auto it = ArcsinChCyTHP.begin(); it != ArcsinChCyTHP.end(); it++){
            AllArcs.push_back(make_pair(it->first, it->second.size()));
        }
        //Sort elements
        sort(AllArcs.begin(), AllArcs.end(), sortCovElms);
        for (int i = 0; i < AllArcs.size(); i++){
            if (UsedArcs < MaxArcFailures){
                for (int j = 0; j < Arcsin1stStageSol.size();j++){//If AllArcs[i] is in Arcsin1stStageSol
                    if (Arcsin1stStageSol[j] == AllArcs[i].first){
                        UsedArcs++;
                        scenarioHeuristics.push_back(AllArcs[i].first);
                        break;
                    }
                }
            }
        }
    }

    //NodeCyChTPH
    for (auto it = NodeCyChTPH.begin(); it != NodeCyChTPH.end(); it++){
        AllVex.push_back(make_pair(it->first, it->second.size()));
    }
    //Sort elements
    sort(AllVex.begin(), AllVex.end(), sortPairs);
    for (int i = 0; i < AllVex.size(); i++){
        if (UsedVertices < MaxVertexFailures){
            for (int j = 0; j < ListSelVex.size();j++){//If AllArcs[i] is in ListSelVex
                if (ListSelVex[j] == AllVex[i].first){
                    UsedVertices++;
                    scenarioHeuristics.push_back(make_pair(-1, AllVex[i].first));
                    break;
                }
            }
        }
    }
}
bool IsxinChain(int v, Chain c){
    for (int i = 0; i < c.Vnodes.size(); i++){
        if (c.Vnodes[i].vertex == v){
            return true;
        }
    }
    return false;
}
bool Problem::UnMVtxdueToVtx(vector<int>& FailedVertices, vector<pair<int,int>>& FailedArcs, vector<int> vinFirstStage, pair<int,int> origin){
    //Build Adjacency List
    IloNumArray2 AdjaList (env);
    //Check whether vertex it->first fails naturally due to another vertex failure
    //Create modified Adjacency List
    AdjaList = BuildAdjaListVtxCycles(FailedVertices, FailedArcs, vinFirstStage,
                                      origin);
    //Call SubCycleFinder
    vector<Cycles> ListC;
    if (origin.first != -1){
        ListC = SubCycleFinder(env, AdjaList, origin.first);
    }
    else{
        ListC = SubCycleFinder(env, AdjaList, origin.second);
    }
    //If no cycle then:
    if (ListC.size() == 0){
        //Build AdjacencyList to find chains
        AdjaList = BuildAdjaListVtxChains(FailedVertices, FailedArcs, vinFirstStage, origin);
        
        //Check whether there's a feasible path from an NDD to vertex i
        vector<int>ChainStarters;
//        if (origin.first != -1){
//            ChainStarters.push_back(origin.first);
//        }
//        else{
        ChainStarters.push_back(origin.second);
//        }
        vector<vChain> VertexinSolChain;
        vector<int>Altruists;
        vector<int>Vertices;
        if (RecoursePolicy == "Among" || RecoursePolicy == "BackArcs"){
            for (int i = 0; i < vinFirstStage.size(); i++){
                if (vinFirstStage[i] > Pairs - 1) Altruists.push_back(i);
            }
        }
        else{//Full
            for (int i = 0; i < Nodes; i++) {
                Vertices.push_back(i);
                if (i > Pairs - 1){
                    Altruists.push_back(i);
                }
            }
        }
        if (RecoursePolicy == "Full"){
            InitializeVertexinSolChain(Vertices, VertexinSolChain, AdjaList);
        }
        else{
            InitializeVertexinSolChain(vinFirstStage, VertexinSolChain, AdjaList);
        }
        vector<Chain>ChainResult;
        ChainResult = FindChains(VertexinSolChain, Altruists, ChainStarters, true);
        if (ChainResult[0].Vnodes.size() == 0){
            return false;
        }
        else{
            return true;
        }
    }
    else{
        return true;
    }
            
}
IloNumArray2 Problem::BuildAdjaListVtxCycles(vector<int> delete_vertex, vector<pair<int, int>> delete_arc, vector<int> vinFirstStage, pair<int,int>origin){
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    
    if (RecoursePolicy == "Full"){
        for (int i = 0; i < AdjacencyList.getSize(); i++){
            AdjaList[i] = IloNumArray(env);
            if (i == origin.first){
                AdjaList[i].add(origin.second + 1);
            }
            else{
                AdjaList[i] = AdjacencyList[i];
            }
        }
    }
    else if (RecoursePolicy == "Among" || RecoursePolicy == "BackArcs"){
        for (int i = 0; i < vinFirstStage.size(); i++){
            if (vinFirstStage[i] == origin.first){
                AdjaList[vinFirstStage[i]].add(origin.second + 1);
            }
            else{
                AdjaList[vinFirstStage[i]] = AdjacencyList[vinFirstStage[i]];
            }
        }
    }
    
    //For all policies
    if (delete_vertex.size() > 0){
        for(int i = 0; i < delete_vertex.size(); i++) AdjaList[delete_vertex[i]] = IloNumArray(env);
    }
    if (delete_arc.size() > 0){
        for (int i = 0; i < delete_arc.size(); i++){
            int val = delete_arc[i].second + 1;
            IloNumArray newrow (env);
            for (int j = 0; j < AdjaList[delete_arc[i].first].getSize(); j++){
                if (AdjaList[delete_arc[i].first][j] != val){
                    newrow.add(AdjaList[delete_arc[i].first][j]);
                }
            }
            //Replace new row into AdjaList
            AdjaList[delete_arc[i].first] = newrow;
        }
    }
    
    return AdjaList;
}
bool isArcTobeDeleted(vector<pair<int, int>> delete_arc, pair<int, int> arc){
    for (int i = 0; i < delete_arc.size(); i++){
        if (arc == delete_arc[i])  return true;
    }
    return false;
}
IloNumArray2 Problem::BuildAdjaListVtxChains(vector<int> delete_vertex, vector<pair<int, int>> delete_arc, vector<int> vinFirstStage, pair<int,int> origin){
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    
    if (RecoursePolicy == "Full"){
        for (int i = 0; i < PredMap.size(); i++){
            AdjaList[i] = IloNumArray(env);
            if ((i == origin.second) && (origin.first != -1)){
                AdjaList[i].add(origin.first + 1);
            }
            else{
                if (IsxinStack(i, delete_vertex) == false){
                    for (int j = 0; j < PredMap[i].size(); j++){
                        if (isArcTobeDeleted(delete_arc, make_pair(PredMap[i][j].first, i)) == false){
                            AdjaList[i].add(PredMap[i][j].first + 1);
                        }
                    }
                }
            }
        }
    }
    else if (RecoursePolicy == "Among" || RecoursePolicy == "BackArcs"){
        for (int i = 0; i < PredMap.size(); i++){
            AdjaList[i] = IloNumArray(env);
            if ((i == origin.second) && (origin.first != -1)){
                AdjaList[i].add(origin.first + 1);
            }
            else{
                if (IsxinStack(i, delete_vertex) == false && IsxinStack (i, vinFirstStage) == true){
                    for (int j = 0; j < PredMap[i].size(); j++){
                        //int v = AdjacencyList[PredMap[i][j].first][PredMap[i][j].second];
                        if (isArcTobeDeleted(delete_arc, make_pair(PredMap[i][j].first,i)) == false && IsxinStack (PredMap[i][j].first, vinFirstStage) == true){
                            AdjaList[i].add(PredMap[i][j].first + 1);
                        }
                    }
                }
            }
        }
    }
    
    return AdjaList;
}
////////////////////Column Generation//////////////////
bool Problem::ColumnGeneration(map<int,bool>&ub_tcyvar, map<int,bool>&ub_tchvar){
    
    //Create model
    IloEnv env;
    IloModel ColGen(env);
    IloCplex Colcplex(ColGen);
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    if (LeftTime < 0){
        Print2ndStage("TimeOut");
    }
    Colcplex.setOut(env.getNullStream());
    Colcplex.setParam(IloCplex::RootAlg, 1);
    Colcplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
    Colcplex.setParam(IloCplex::Param::TimeLimit, LeftTime);
    Colcplex.setParam(IloCplex::Param::Threads, 1);
    clock_t tStartMP = clock();
    tStartMP = clock();
    
    //Create variables
    IloNumVar y(env, 0, 1, ILOFLOAT);
    IloNumVarArray zcol(env, 0, 0, 1, ILOFLOAT);//AllCycles.size()
    
    //Create variables
    int pimany = 0;
    if (ChainLength == 0){
        pimany = int(Pairs);
    }
    else{
        pimany = int(Nodes);
    }
    
    //Arrays for constraints and objective
    IloRangeArray onecycle(env);
    IloObjective Obj(env);
    //Objective
    Obj = IloAdd(ColGen, IloMaximize(env,-10000*y));
    for (int i = 0; i < pimany; i++){
        string name = "c." + to_string(i);
        const char* cName = name.c_str();
        onecycle.add(IloRange(env, -IloInfinity, y,1, cName));
    }
    ColGen.add(onecycle);
    
    double UpperBound = 0;
    IloBool NewCycleAdded = true;
    IloNumArray solpi (env,pimany);
    vector<int>ColumnsAdded;
    vector<pair<int,double>>dualOrder;
    map<int,int>cycles;//zcol variable, cycle in Cycles2ndStage
    map<int,int>chains;//zcol variable, chain in Chains2ndStage
    vector<int>whichChains;
    int half;
    if (ChainLength <= 3){
        for (int i = 0; i < Chains2ndStage.size(); i++) whichChains.push_back(i);
    }
    else{
        half = Chains2ndStage.size()/2;
        for (int i = 0; i < half; i++){
            whichChains.push_back(i);
        }
    }
    int manyCH = 10;
    if (ChainLength >= 4) manyCH = 5;
    
    while(NewCycleAdded == true){

        Colcplex.setOut(env.getNullStream());
        clock_t tStartMP2 = clock();
        tStartMP2 = clock();
        Colcplex.solve();
        
        if (Colcplex.getStatus() == IloAlgorithm::Infeasible) {
            cout << "Infeasible";
        }
        else{
            //cplex.exportModel("LagProb.lp");
            //Get objective value
            UpperBound = Colcplex.getObjValue();
            //cout << "Upper Bound:" << UpperBound << endl;
            
            //Get dual variables
            Colcplex.getDuals(solpi,onecycle);
//            dualOrder.clear();
//            for (int i = 0; i < solpi.getSize(); i++){
//                dualOrder.push_back(make_pair(i, solpi[i]));
//            }
//            //Sort duals
//            sort(dualOrder.begin(), dualOrder.end(), sortduals);
            
            int naddedCH = 0, naddedCY = 0;
                for (int k = 0; k < whichChains.size(); k++){
                    auto it = ub_tchvar.find(k);
                    double dualWeight = 0, Total = 0;
                    if (it == ub_tchvar.end()){
                        if (THP_Bound == "Strong" || (THP_Method == "DoubleCovering" && THP_Bound == "NoBound")){
                            dualWeight = (Chains2ndStage[whichChains[k]].AccumWeight*Nodes + 1);
                        }
                        else if (THP_Bound == "Simple" || (THP_Method == "Covering" && THP_Bound == "NoBound")){
                            dualWeight = Chains2ndStage[whichChains[k]].AccumWeight;
                        }
                    }
                    else{
                        if (THP_Bound == "Strong" || (THP_Method == "DoubleCovering" && THP_Bound == "NoBound")){
                            dualWeight = 1;
                        }
                    }
                    Total = dualWeight;
                    for (int j = 0; j < Chains2ndStage[whichChains[k]].Vnodes.size(); j++){
                        Total-= solpi[Chains2ndStage[whichChains[k]].Vnodes[j].vertex];
                    }
                    if (Total > 0.005){
                        naddedCH++;
                        //Create chain column
                        IloNumColumn col(env);
                        //Update constraints
                        for (int j = 0; j < Chains2ndStage[whichChains[k]].Vnodes.size(); j++){
                            col+= onecycle[Chains2ndStage[whichChains[k]].Vnodes[j].vertex](1);
                        }
                        //Add weight of new column
                        col += Obj(dualWeight);
                        //Create new chain variable
                        zcol.add(IloNumVar(col, 0, 1));
                        string name = "zcol." + to_string(zcol.getSize());
                        zcol[zcol.getSize() - 1].setName(name.c_str());
                        //Store new variable
                        chains[int(zcol.getSize()) - 1] = whichChains[k];
                        col.end();
                        if (naddedCH >= manyCH) break;
                    }
                    
                }
                for (int k = 0; k < Cycles2ndStage.size(); k++){
                    auto it = ub_tcyvar.find(k);
                    double dualWeight = 0, Total = 0;
                    if (it == ub_tcyvar.end()){
                        if (THP_Bound == "Strong" || (THP_Method == "DoubleCovering" && THP_Bound == "NoBound")){
                            dualWeight = (Cycles2ndStage[k].get_Many()*Nodes + 1);
                        }
                        else if (THP_Bound == "Simple" || (THP_Method == "Covering" && THP_Bound == "NoBound")){
                            dualWeight = Cycles2ndStage[k].get_Many();
                        }
                    }
                    else{
                        if (THP_Bound == "Strong" || (THP_Method == "DoubleCovering" && THP_Bound == "NoBound")){
                            dualWeight = 1;
                        }
                    }
                    Total = dualWeight;
                    for (int j = 0; j < Cycles2ndStage[k].get_c().size(); j++){
                        Total-= solpi[Cycles2ndStage[k].get_c()[j]];
                    }
                    if (Total > 0.005){
                        naddedCY++;
                        //Create chain column
                        IloNumColumn col(env);
                        //Update constraints
                        for (int j = 0; j < Cycles2ndStage[k].get_c().size(); j++){
                            col+= onecycle[Cycles2ndStage[k].get_c()[j]](1);
                        }
                        //Add weight of new column
                        col += Obj(dualWeight);
                        //Create new chain variable
                        zcol.add(IloNumVar(col, 0, 1));
                        string name = "zcol." + to_string(zcol.getSize());
                        zcol[zcol.getSize() - 1].setName(name.c_str());
                        //Store new variable
                        cycles[int(zcol.getSize()) - 1] = k;
                        col.end();
                        if (naddedCY >= 10) break;
                    }
                }
            //}

            if (naddedCY == 0 && naddedCH == 0){
                NewCycleAdded = false;
                if (whichChains.size() !=  Chains2ndStage.size()){
                    NewCycleAdded = true;
                    for (int i = half; i < Chains2ndStage.size(); i++){
                            whichChains.push_back(i);
                    }
                }
            }
        }
    }
    UpperBound = floor (UpperBound + 0.1);
    //Get feasible solution
    Colcplex.setParam(IloCplex::Param::TimeLimit, 100);
    //Transform variables into integers
    ColGen.add(IloConversion(env, zcol, ILOBOOL));
    //Set start MIP solution
    Colcplex.solve();
    double objCol = Colcplex.getObjValue();
    if (objCol == UpperBound){
        //Retrieve solution
        TPMIP_Obj = UpperBound;
        IloNumArray zsol(env, zcol.getSize());
        tcysolColGen = IloNumArray(env, Cycles2ndStage.size());
        tchsolColGen = IloNumArray(env, Chains2ndStage.size());
        Colcplex.getValues(zsol,zcol);
        for (int i = 0; i < zsol.getSize(); i++){
            if (zsol[i] > 0.9){
                auto it = cycles.find(i);
//                cout << endl;
                if (it!= cycles.end()){
                    tcysolColGen[cycles[i]] = 1;
//                    for (int j = 0; j < Cycles2ndStage[cycles[i]].get_c().size();j++){
//                        cout << Cycles2ndStage[cycles[i]].get_c()[j] <<"\t";
//                    }
                }
                else{
                    tchsolColGen[chains[i]] = 1;
//                    for (int j = 0; j < Chains2ndStage[chains[i]].Vnodes.size();j++){
//                        cout << Chains2ndStage[chains[i]].Vnodes[j].vertex << "\t";
//                    }
                }
            }
        }
        
        ColGen.end();
        Colcplex.end();
        return true;
    }
    
    ColGen.end();
    Colcplex.end();
    return false;
}
