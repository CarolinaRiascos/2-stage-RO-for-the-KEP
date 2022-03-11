//
//  ThirdPhase.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "ThirdPhase.hpp"
void Problem::THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<int>&vinFirstStage){
    //Get scenario
    
    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
    //Create model
    mTHPMIP = IloModel (env);
    cplexmTHPMIP = IloCplex(mTHPMIP);
    cplexmTHPMIP.setParam(IloCplex::Param::Threads, 1);
    cplexmTHPMIP.setOut(env.getNullStream());
    
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
    
    Ite2ndS = 0;
    while(true){
        tEnd2ndS = (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
        if (tEnd2ndS > 1800){
            status = "Incomplete";
            break;
        }
        Ite2ndS++;
        //Solve formulation
        //cplexmTHPMIP.exportModel("mTHPMIP.lp");
        tStartReco = clock();
        cplexmTHPMIP.solve();
        tTotalReco += (clock() - tStartReco)/double(CLOCKS_PER_SEC);
        if (cplexmTHPMIP.getStatus() == IloAlgorithm::Infeasible){
            cout << "S.O.S. This should not happen." << endl;
        }
        else{
            TPMIP_Obj = cplexmTHPMIP.getObjValue();
                        
            IloNumArray tcysol(env);
            IloNumArray tchsol(env);
            cplexmTHPMIP.getValues(tcysol,tcyvar);
            cplexmTHPMIP.getValues(tchsol,tchvar);

            Get3rdStageSol(Cycles3rdSol, Chains3rdSol, tcysol, tchsol);
            
            if (THP_Method == "Covering"){
                if (ThisWork(tcysol, tchsol, vinFirstStage) == true) break;
            }else{//Literature's approach
                if (Literature(tcysol, tchsol) == true) break;
            }
            
            tcysol.end();
            tchsol.end();
        }
    }
    if (tEnd2ndS <= 1800) status = "Complete";
    Print2ndStage();
    
    //Check whether a new scenario was found
    if (SPMIP_Obj < FPMIP_Obj){
        //Send scenario to first phase
        //Check for FailedArcs and FailedVertices when back to the first stage
        
    }else if (SPMIP_Obj == FPMIP_Obj) {
        //Robust solution found!
    }
    else{
        cout << "S.O.S. This should never happen";
    }
    
}
bool Problem::ThisWork(IloNumArray& tcysol, IloNumArray& tchsol, vector<int>&vinFirstStage){
    
    if (TPMIP_Obj < LOWEST_TPMIP_Obj){
        LOWEST_TPMIP_Obj = TPMIP_Obj;
        OptFailedArcs = FailedArcs;
        OptFailedVertices = FailedVertices;
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
    GetAtLeastOneFails(Cycles3rdSol, Chains3rdSol, vinFirstStage);
    
    //Resolve 2nd. Phase
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    bool runH = false;//Heuristcs2ndPH();
    if (runH == false){
        cplexGrandSubP.solve();

        if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
            cout << "2nd. stage solved" << endl;
            SPMIP_Obj = LOWEST_TPMIP_Obj;
            return true;
        }

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
                    if (AdjacencyList[scenarioHeuristics[j].first][i] == scenarioHeuristics[j].second){
                        arc_sol[scenarioHeuristics[j].first][i] = 1;
                        break;
                    }
                }
                
            }
        }
    }
    GetScenario(arc_sol, vertex_sol); //Update FailedArcs and FailedVertices
    
    //Change variables' UB in 3rd phase
    map<int,bool>ub_tcyvar; // <cycle number, 0 or 1>
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
    
    map<int,bool>ub_tchvar; // <chain number, 0 or 1>
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
    
    map<int,bool>ub_tcyvar; // <cycle number, 0 or 1>
    map<int,bool>ub_tchvar; // <cycle number, 0 or 1>
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
void Problem::GetAtLeastOneFails(vector<Cycles>&Cycles3rdSol, vector<Chain>&Chains3rdSol, vector<int>&vinFirstStage){
    IloExpr expr(env);
    
    
    for (int i = 0; i < Cycles3rdSol.size(); i++){
        for (int j = 0; j < Cycles3rdSol[i].get_c().size(); j++){
            int u = Cycles3rdSol[i].get_c()[j];
            auto it = Elms2ndPhase.find(make_pair(-1,u));
            if (it == Elms2ndPhase.end()){
                Elms2ndPhase[make_pair(-1,u)] = coveringElements();
            }
            Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
            expr+= vertex[u];
            if (j == Cycles3rdSol[i].get_c().size() - 1){
                int s = mapArcs[make_pair(u,Cycles3rdSol[i].get_c()[0])];
                expr+= arc[u][s];
                auto it = Elms2ndPhase.find(make_pair(u,Cycles3rdSol[i].get_c()[0]));
                if (it == Elms2ndPhase.end()){
                    Elms2ndPhase[make_pair(u,Cycles3rdSol[i].get_c()[0])] = coveringElements();
                }
                Elms2ndPhase[make_pair(u,Cycles3rdSol[i].get_c()[0])].add_const(int(AtLeastOneFails.getSize()));
            }
            else{
                int s = mapArcs[make_pair(u,Cycles3rdSol[i].get_c()[j + 1])];
                expr+= arc[u][s];
                auto it = Elms2ndPhase.find(make_pair(u,Cycles3rdSol[i].get_c()[j + 1]));
                if (it == Elms2ndPhase.end()){
                    Elms2ndPhase[make_pair(u,Cycles3rdSol[i].get_c()[j + 1])] = coveringElements();
                }
                Elms2ndPhase[make_pair(u,Cycles3rdSol[i].get_c()[j + 1])].add_const(int(AtLeastOneFails.getSize()));
            }
        }
    }
    for (int i = 0; i < Chains3rdSol.size(); i++){
        for (int j = 0; j < Chains3rdSol[i].Vnodes.size(); j++){
            int u = Chains3rdSol[i].Vnodes[j].vertex;
            auto it = Elms2ndPhase.find(make_pair(-1,u));
            if (it == Elms2ndPhase.end()){
                Elms2ndPhase[make_pair(-1,u)] = coveringElements();
            }
            Elms2ndPhase[make_pair(-1,u)].add_const(int(AtLeastOneFails.getSize()));
            expr+= vertex[u];
            if (j <= Chains3rdSol[i].Vnodes.size() - 2){
                int s = mapArcs[make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex)];
                expr+= arc[u][s];
                auto it = Elms2ndPhase.find(make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex));
                if (make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex) == make_pair(34, 34)){
                    cout << "hi";
                }
                if (it == Elms2ndPhase.end()){
                    Elms2ndPhase[make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex)] = coveringElements();
                }
                Elms2ndPhase[make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex)].add_const(int(AtLeastOneFails.getSize()));
            }
        }
    }
    
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
        sort(RecoSolCovering.back().begin(), RecoSolCovering.back().end(), sortint);
    }

    if (RecoTotalWCovering.back() > LOWEST_TPMIP_Obj){
        RHS = Update_RHS_Covering(int (RecoSolCovering.size() - 1));
    }
    //Add new constraint to Const2ndPhase
    Const2ndPhase.push_back(RHS);
    
    string name = "AtLeastOneFails_" + to_string(AtLeastOneFails.getSize() + 1);
    const char* cName = name.c_str();
    //cout << expr << endl;
    AtLeastOneFails.add(IloRange(env, RHS, expr, IloInfinity, cName));
    GrandSubProb.add(AtLeastOneFails[AtLeastOneFails.getSize() - 1]);
    expr.end();
}
bool sortint(double& c1, double& c2){
    return (c1 > c2);
}
bool sortCovElms(pair< pair<int,int>, int>& c1, pair< pair<int,int>, int>& c2){
    return (c1.second > c2.second);
}
int Problem::Update_RHS_Covering(int row){
    int i = 0, RHS = 0, accum = 0;
        
    while (true){
        if (RecoTotalWCovering[row] - RecoSolCovering[row][i] - accum > LOWEST_TPMIP_Obj){
            accum+= RecoSolCovering[row][i];
            i++;
        }
        else{
            RHS = i + 1;
            break;
        }
    }
    if (i >= 2) VI_I++;
    
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
    if (THP_Method != "Literature"){
        for(int i = 0; i < Cycles2ndStage.size(); i++){
            int n = int(Cycles2ndStage[i].get_Many());
                expr += n*tcyvar[i];
        }
        for(int i = 0; i < Chains2ndStage.size(); i++){
            int n = int(Chains2ndStage[i].AccumWeight);
                expr += n*tchvar[i];
        }
        return expr;
    }else{
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
    if (THP_Method != "Literature"){
        ub_tcyvar = GetUB_tcyvar(FailedArcs, FailedVertices);
    }
    
    IloNumVarArray var(env, Cycles2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        if (THP_Method != "Literature"){
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
    if (THP_Method != "Literature"){
        ub_tchvar = GetUB_tchvar(FailedArcs, FailedVertices);
    }
    
    IloNumVarArray var(env, Chains2ndStage.size());
    for (int i = 0; i < var.getSize(); i++){
        if (THP_Method != "Literature"){
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
bool Problem::Heuristcs2ndPH(){
    bool another = false, keepgoing = true;
    int UsedArcs = 0;
    int UsedVertices = 0, ele, ite = 0;
    vector<pair<int,int>>tabuList;
    
    while(ite <= 20 && keepgoing == true){
        scenarioHeuristics.clear();
        tabuList.clear();
        for (auto it = Elms2ndPhase.begin(); it != Elms2ndPhase.end(); it++) it->second.set_state(false);
        for (int i = 0; i < Const2ndPhase.size(); i++){
            Const2ndPhase[i].deleteEls();
        }
        UsedArcs = 0;
        UsedVertices = 0;
        ite++;
        while (true){
            vector<pair< pair<int,int>, int>>auxCov;
            for (auto it = Elms2ndPhase.begin(); it != Elms2ndPhase.end(); it++){
                if (it->first.first == -1 && UsedVertices < MaxVertexFailures && Elms2ndPhase[it->first].get_state() == false){
                    auxCov.push_back(make_pair(it->first, it->second.get_coversize()));
                }
                if (it->first.first != -1 && UsedArcs < MaxArcFailures && Elms2ndPhase[it->first].get_state() == false){
                    auxCov.push_back(make_pair(it->first, it->second.get_coversize()));
                }
            }
            //If no Elements to pick from
            if (auxCov.size() == 0){
                break;
            }
            //Sort elements
            sort(auxCov.begin(), auxCov.end(), sortCovElms);
            while(true){
                int search = auxCov[0].second;
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
            // Add element to scenario
            scenarioHeuristics.push_back(auxCov[ele].first);
            //Seize the element
            Elms2ndPhase[auxCov[ele].first].set_state(true);
            if (auxCov[ele].first.first == -1){//It is a vertex
                UsedVertices++;
            }
            else{
                UsedArcs++;
            }
            //Cover constraints
            for (int i = 0; i < Elms2ndPhase[auxCov[ele].first].get_coveredconsts().size(); i++){
                int e = Elms2ndPhase[auxCov[ele].first].get_coveredconsts()[i];
                Const2ndPhase[e].add_cover(auxCov[ele].first);
            }
            keepgoing = false;
            //Check whether there is an unsatisfied constraint
            for (int i = 0; i < Const2ndPhase.size(); i++){
                if (Const2ndPhase[i].get_coversize() < Const2ndPhase[i].get_RHS()){
                    keepgoing = true;
                }
            }
            if (keepgoing == false) break;
        }
    }
    if (keepgoing == true){
        return false;
    }
    else{
        return true;
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
bool Problem::UnMVtxdueToVtx(map<int, bool>& FailedVertices, vector<int> vinFirstStage){
    //Build Adjacency List
    IloNumArray2 AdjaList (env);
    bool thereischain = false;
    //Check whether vertex it->first fails naturally due to another vertex failure
    for (auto it = FailedVertices.begin(); it!= FailedVertices.end(); it++){
        vector<int>others;
        for (auto it2 = FailedVertices.begin(); it2 != FailedVertices.end(); it2++){
            if (it->first != it2->first){
                others.push_back(it2->first);
            }
        }
        //Create modified Adjacency List
        AdjaList = BuildAdjaListVtxCycles(others, vector<pair<int,int>>(), vinFirstStage);
        //Call SubCycleFinder
        int origin = it->first;
        vector<Cycles> ListC;
        ListC = SubCycleFinder(env, AdjaList, origin);
        //If no cycle then:
        if (ListC.size() == 0){
            //Build AdjacencyList to find chains
            AdjaList = BuildAdjaListVtxChains(others, vector<pair<int,int>>(), vinFirstStage);
            
            //Check whether there's a feasible path from an NDD to vertex i
            vector<int>ChainStarters; ChainStarters.push_back(it->first);
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
            if (ChainResult.size() == 0){
                cout << "Hi";
            }
        }
            
    }
    return thereischain;
}
IloNumArray2 Problem::BuildAdjaListVtxCycles(vector<int> delete_vertex, vector<pair<int, int>> delete_arc, vector<int> vinFirstStage){
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    
    if (RecoursePolicy == "Full"){
        for (int i = 0; i < AdjacencyList.getSize(); i++){
            AdjaList[i] = AdjacencyList[i];
        }
    }
    else if (RecoursePolicy == "Among" || RecoursePolicy == "BackArcs"){
        for (int i = 0; i < vinFirstStage.size(); i++){
            AdjaList[vinFirstStage[i]] = AdjacencyList[vinFirstStage[i]];
        }
    }
    
    //For all policies
    if (delete_vertex.size() > 0){
        for(int i = 0; i < delete_vertex.size(); i++) AdjaList[delete_vertex[i]] = IloNumArray(env);
    }
    if (delete_arc.size() > 0){
        for (int i = 0; i < delete_arc.size(); i++){
            int val = AdjacencyList[delete_arc[i].first][delete_arc[i].second];
            IloNumArray newrow (env);
            for (int j = 0; j < AdjaList[delete_arc[i].first].getSize(); j++){
                if (AdjaList[delete_arc[i].first][j] != val){
                    newrow.add(AdjaList[i]);
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
IloNumArray2 Problem::BuildAdjaListVtxChains(vector<int> delete_vertex, vector<pair<int, int>> delete_arc, vector<int> vinFirstStage){
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    
    if (RecoursePolicy == "Full"){
        for (int i = 0; i < PredMap.size(); i++){
            AdjaList[i] = IloNumArray(env);
            if (IsxinStack(i, delete_vertex) == false){
                for (int j = 0; j < PredMap[i].size(); j++){
                    if (isArcTobeDeleted(delete_arc, make_pair(PredMap[i][j].first, i)) == false){
                        AdjaList[i].add(PredMap[i][j].first);
                    }
                }
            }
        }
    }
    else if (RecoursePolicy == "Among" || RecoursePolicy == "BackArcs"){
        for (int i = 0; i < PredMap.size(); i++){
            if (IsxinStack(i, delete_vertex) == false && IsxinStack (i, vinFirstStage) == true){
                for (int j = 0; j < PredMap[i].size(); j++){
                    int v = AdjacencyList[PredMap[i][j].first][PredMap[i][j].second];
                    if (isArcTobeDeleted(delete_arc, make_pair(v,i)) == false && IsxinStack (v, vinFirstStage) == true){
                        AdjaList[i].add(v);
                    }
                }
            }
        }
    }
    
    return AdjaList;
}
