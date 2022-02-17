//
//  ThirdPhase.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "ThirdPhase.hpp"
void Problem::THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<int>&ListSelVertices){
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
    bool RemoveExcess = false;
    while(true){
        Ite2ndS++;
        //Solve formulation
        //cplexmTHPMIP.exportModel("mTHPMIP.lp");
        cplexmTHPMIP.solve();
        if (cplexmTHPMIP.getStatus() == IloAlgorithm::Infeasible){
            cout << "S.O.S. This should not happen." << endl;
        }
        else{
            TPMIP_Obj = cplexmTHPMIP.getObjValue();
            if (TPMIP_Obj < LOWEST_TPMIP_Obj){
                LOWEST_TPMIP_Obj = TPMIP_Obj;
                OptFailedArcs = FailedArcs;
                OptFailedVertices = FailedVertices;
            }
            IloNumArray tcysol(env);
            IloNumArray tchsol(env);
            cplexmTHPMIP.getValues(tcysol,tcyvar);
            cplexmTHPMIP.getValues(tchsol,tchvar);
            Get3rdStageSol(Cycles3rdSol, Chains3rdSol, tcysol, tchsol);
            if (THP_Method == "Covering" || THP_Method == "BoundedCovering"){
                if (ThisWork(tcysol, tchsol, ListSelVertices, RemoveExcess) == true) break;
            }else{//Literature's approach
                if (Literature(tcysol, tchsol) == true) break;
            }
            
            tcysol.end();
            tchsol.end();
        }
    }
    tEnd2ndS = (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
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
bool Problem::ThisWork(IloNumArray& tcysol, IloNumArray& tchsol, vector<int>&ListSelVertices, bool RemoveExcess){
    

    if (THP_Method == "BoundedCovering"){
        //Check whether the cycles are contained in Cycles2ndTo3rd and Chains2ndTo3rd
        AddNewCols3rdTo2nd (tcysol, tchsol, Cycles2ndTo3rd, Chains2ndTo3rd, Cycles3rdTo2nd, Chains3rdTo2nd, Cyclenewrow2ndPH, Chainnewrow2ndPH, tcysol3rd, tchsol3rd, Cycles2ndStage, Chains2ndStage);

        if (tcysol3rd.size() != 0 || tchsol3rd.size() != 0){
            //Add cols and constraints to SPMIP
            for (int i = 0; i < tcysol3rd.size(); i++){
                
                //OldSelVert2ndPH_CY
                vector<int> OldSelVert2ndPH;
                //OldSelVert2ndPH = ModifyOldSelectVex_CY(tcysol3rd[i], ListSelVertices, Cycles2ndStage);
                //OldActiveCCSubP_CY
                vector<int> OldActiveCY;
                OldActiveCY = ModifyOldActiveCCSubP_CY(tcysol3rd[i], Cycles2ndTo3rd, Cyclenewrow2ndPH, Cycles2ndStage);
                //OldActiveCCSubP_CH
                vector<int> OldActiveCYtoCH;
                OldActiveCYtoCH = ModifyOldActiveCCSubP_CYtoCH(tcysol3rd[i], Chains2ndTo3rd, Chainnewrow2ndPH, Cycles2ndStage);

                //Create new cycle column
                IloNumColumn col(env);
                
                //Update Beta constraint
                col+= ConsBeta[0](-Cycles2ndStage[tcysol3rd[i]].get_Many());
                //UpdateSelVert2ndPH
//                for (int j = 0; j < OldSelVert2ndPH.size(); j++){
//                    col+= SelVert2ndPH[OldSelVert2ndPH[j]](1);
//                }
                //UpdateActiveCCSubP_CY
                for (int j = 0; j < OldActiveCY.size(); j++){
                    col+= ActiveCCSubP_CY[OldActiveCY[j]](1);
                }
                //UpdateActiveCCSubP_CH
                for (int j = 0; j < OldActiveCYtoCH.size(); j++){
                    col+= ActiveCCSubP_CH[OldActiveCYtoCH[j]](1);
                }
                //Update selOnce
                for (int j = 0; j < Cycles2ndStage[tcysol3rd[i]].get_c().size(); j++){
                    auto it = CycleNodeTPH.find(Cycles2ndStage[tcysol3rd[i]].get_c()[j]);
                    if (it != CycleNodeTPH.end()){
                        col+= Once[Cycles2ndStage[tcysol3rd[i]].get_c()[j]](1);
                    }
                }
                
                //Create new variable
                cyvar.add(IloIntVar(col, 0, 1));
                string name = "x." + to_string(Cyclenewrow2ndPH + i + 1);
                const char* varName = name.c_str();
                cyvar[cyvar.getSize() - 1].setName(varName);
                cyvar[cyvar.getSize() - 1].setBounds(0, 1);
                //NewActiveCCSubP_CY
                GetNewIloRangeCY3rd(tcysol3rd[i], Cycles2ndStage);
                //IfNoFailuresCY
                IfNoFailuresCY(Cycles2ndTo3rd, Cycles3rdTo2nd[tcysol3rd[i]]);
            }
            for (int i = 0; i < tchsol3rd.size(); i++){
                //OldSelVert2ndPH
                //vector<int> OldSelVert2ndPH;
                //OldSelVert2ndPH = ModifyOldSelectVex_CH(tchsol3rd[i], ListSelVertices, Chains2ndStage);
                //OldActiveCCSubP_CH
                vector<int> OldActiveCH;
                OldActiveCH = ModifyOldActiveCCSubP_CH(tchsol3rd[i], Chains2ndTo3rd, Chainnewrow2ndPH, Chains2ndStage);
                //OldActiveCCSubP_CH
                vector<int> OldSelVert2ndPH_CHtoCY;
                //OldActiveCCSubP_CHtoCY
                OldSelVert2ndPH_CHtoCY = ModifyOldActiveCCSubP_CHtoCY(tchsol3rd[i], Cycles2ndTo3rd, Cyclenewrow2ndPH, Chains2ndStage);
                
                //Create new chain column
                IloNumColumn col(env);
                
                //Update Beta constraint
                col+= ConsBeta[0](-Chains2ndStage[tchsol3rd[i]].AccumWeight);
                //UpdateSelVert2ndPH
//                for (int j = 0; j < OldSelVert2ndPH.size(); j++){
//                    col+= SelVert2ndPH[OldSelVert2ndPH[j]](1);
//                }
                //UpdateActiveCCSubP_CH
                for (int j = 0; j < OldActiveCH.size(); j++){
                    col+= ActiveCCSubP_CH[OldActiveCH[j]](1);
                }
                //UpdateActiveCCSubP_CY
                for (int j = 0; j < OldSelVert2ndPH_CHtoCY.size(); j++){
                    col+= ActiveCCSubP_CY[OldSelVert2ndPH_CHtoCY[j]](1);
                }
                
                //Update selOnce
                for (int j = 0; j < Chains2ndStage[tchsol3rd[i]].Vnodes.size(); j++){
                    auto it = ChainNodeTPH.find(Chains2ndStage[tchsol3rd[i]].Vnodes[j].vertex);
                    if (it != ChainNodeTPH.end()){
                        col+= Once[Chains2ndStage[tchsol3rd[i]].Vnodes[j].vertex](1);
                    }
                }
                
                //Create new variable
                chvar.add(IloIntVar(col, 0, 1));
                string name = "y." + to_string(Chainnewrow2ndPH + i + 1);
                const char* varName = name.c_str();
                chvar[chvar.getSize() - 1].setName(varName);
                chvar[chvar.getSize() - 1].setBounds(0, 1);
                //NewActiveCCSubP_CH
                GetNewIloRangeCH3rd(tchsol3rd[i], Chains2ndStage);
                //IfNoFailuresCY
                IfNoFailuresCH(Chains2ndTo3rd, Chains3rdTo2nd[tchsol3rd[i]]);
            }
        }
        
        
    }

    if (THP_Method == "Covering"){
        //Add AtLeastOneFails Cut
        GetAtLeastOneFails(Cycles3rdSol, Chains3rdSol);
    }else{
        //Add New Beta Constraint
        GetNewBetaCut(TPMIP_Obj, FailedArcs, FailedVertices);
        GetAtLeastOneFails(Cycles3rdSol, Chains3rdSol);
    }
    
    
    //Resolve 2nd. Phase
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    cplexGrandSubP.solve();
    
    if (THP_Method == "Covering"){
        if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
            cout << "2nd. stage solved" << endl;
            SPMIP_Obj = LOWEST_TPMIP_Obj;
            return true;
        }
    }
    else{
        SPMIP_Obj = cplexGrandSubP.getValue(Beta);
        if (SPMIP_Obj >= TPMIP_Obj){
             cout << "2nd. stage solved" << endl;
            //SPMIP_Obj = LOWEST_TPMIP_Obj;
            //return true;
        }
    }


    if (THP_Method == "BoundedCovering"){
        if (SPMIP_Obj >= TPMIP_Obj){
            cout << "This should not happen" << endl;
        }
        cout << "2nd. Phase Obj: " << SPMIP_Obj << endl;
        
        cyvar_sol2nd = IloNumArray(env, Cycles2ndTo3rd.size());
        chvar_sol2nd = IloNumArray (env, Chains2ndTo3rd.size());
        cplexGrandSubP.getValues(cyvar_sol2nd,cyvar);
        cplexGrandSubP.getValues(chvar_sol2nd,chvar);
        
//        cout << "Sol. 2nd. Phase: " << endl;
//        for (int i = 0; i < cyvar_sol2nd.getSize(); i++){
//            if (cyvar_sol2nd[i] > 0.1) cout << cyvar[i].getName() << endl;
//        }
//        for (int i = 0; i < chvar_sol2nd.getSize(); i++){
//            if (chvar_sol2nd[i] > 0.1) cout << chvar[i].getName() << endl;
//        }
        cyvar_sol2nd.end();
        chvar_sol2nd.end();
    }

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
                Cycles3rdTo2nd[i] = int(Cycles3rdTo2nd.size());
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
                Chains3rdTo2nd[i] = int(Chains3rdTo2nd.size());
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
        SPMIP_Obj = cplexGrandSubP.getValue(Beta);
        if (SPMIP_Obj >= TPMIP_Obj){
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
void Problem::GetAtLeastOneFails(vector<Cycles>&Cycles3rdSol, vector<Chain>&Chains3rdSol){
    map<int, bool>::iterator it1;
    map<pair<int,int>, bool>::iterator it2;
    bool alreadyfails = false;
    IloExpr expr(env);
    
    for (int i = 0; i < Cycles3rdSol.size(); i++){
        for (int j = 0; j < Cycles3rdSol[i].get_c().size(); j++){
            int u = Cycles3rdSol[i].get_c()[j];
            expr+= vertex[u];
            it1 = FailedVertices.find(u);
            if (it1 != FailedVertices.end()) {
                alreadyfails = true;
                break;
            }
            if (j == Cycles3rdSol[i].get_c().size() - 1){
                int s = mapArcs[make_pair(u,Cycles3rdSol[i].get_c()[0])];
                expr+= arc[u][s];
                it2 = FailedArcs.find(make_pair(u,Cycles3rdSol[i].get_c()[0]));
                if (it2 != FailedArcs.end()){
                    alreadyfails = true;
                    break;
                }
            }
            else{
                int s = mapArcs[make_pair(u,Cycles3rdSol[i].get_c()[j + 1])];
                expr+= arc[u][s];
                it2 = FailedArcs.find(make_pair(u,Cycles3rdSol[i].get_c()[j + 1]));
                if (it2 != FailedArcs.end()){
                    alreadyfails = true;
                    break;
                }
            }
        }
        if (alreadyfails == true) break;
    }
    if (alreadyfails == false){
        for (int i = 0; i < Chains3rdSol.size(); i++){
            for (int j = 0; j < Chains3rdSol[i].Vnodes.size(); j++){
                int u = Chains3rdSol[i].Vnodes[j].vertex;
                it1 = FailedVertices.find(u);
                if (it1 != FailedVertices.end()){
                    alreadyfails = true;
                    break;
                }
                expr+= vertex[u];
                if (j <= Chains3rdSol[i].Vnodes.size() - 2){
                    int s = mapArcs[make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex)];
                    expr+= arc[u][s];
                    it2 = FailedArcs.find(make_pair(u,Chains3rdSol[i].Vnodes[j + 1].vertex));
                    if (it2 != FailedArcs.end()){
                        alreadyfails = true;
                        break;
                    }
                }
            }
            if (alreadyfails == true) break;
        }
    }
    if (alreadyfails == false){
        string name = "AtLeastOneFails." + to_string(AtLeastOneFails.getSize() + 1);
        const char* cName = name.c_str();
        //cout << expr << endl;
        AtLeastOneFails.add(IloRange(env, 1, expr, IloInfinity, cName));
        GrandSubProb.add(AtLeastOneFails[AtLeastOneFails.getSize() - 1]);
    }
    expr.end();
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
vector<int>Problem::ModifyOldActiveCCSubP_CHtoCY(int tOnechsol, map<int,int>& Cycles2ndTo3rd, int& Cyclenewrow2ndPH, vector<Chain>&Chains2ndStage){
    vector<int>UpdateActiveCCSubP_CY;
    for (int i = 0; i < Cyclenewrow2ndPH; i++){
        bool found = false;
        for (int j = 0; j < Cycles2ndStage[Cycles2ndTo3rd[i]].get_c().size(); j++){
            int u = Cycles2ndStage[Cycles2ndTo3rd[i]].get_c()[j];
            for (int h = 0; h < Chains2ndStage[tOnechsol].Vnodes.size(); h++){
                int v = Chains2ndStage[tOnechsol].Vnodes[h].vertex;
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
vector<int>Problem::ModifyOldActiveCCSubP_CYtoCH(int tOnecysol, map<int,int>&Chains2ndTo3rd, int&Chainnewrow2ndPH, vector<Cycles>&Cycles2ndStage){
    vector<int>UpdateActiveCCSubP_CH;
    for (int i = 0; i < Chainnewrow2ndPH; i++){
        bool found = false;
        for (int j = 0; j < Chains2ndStage[Chains2ndTo3rd[i]].Vnodes.size(); j++){
            int u = Chains2ndStage[Chains2ndTo3rd[i]].Vnodes[j].vertex;
            for (int h = 0; h < Cycles2ndStage[tOnecysol].get_c().size(); h++){
                int v = Cycles2ndStage[tOnecysol].get_c()[h];
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
                break;
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
                break;
            }
        }
    }
    
    return UpdateOldSelectVex_CH;
}
void Problem::GetNewIloRangeCY3rd(int idx, vector<Cycles>&Cycles2ndStage){
    
    IloExpr ExprArcsVtx (env, 0);
    IloExpr AllccVars (env, 0);
    vector<int>forcycles;
    vector<int>forchains;
    for (int j = 0; j < Cycles2ndStage[idx].get_c().size(); j++){
        int u = Cycles2ndStage[idx].get_c()[j];
        ExprArcsVtx+= vertex[u];
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
            if (CycleNodeSPH[u][k] != Cycles3rdTo2nd[idx]) forcycles.push_back(CycleNodeSPH[u][k]);
            //check whether CycleNodeSPH[Cycles2ndStage[idx].get_c()[j]][k] exixts in Cycles2ndTo3rd
        }
        //cout << Cycles2ndStage[idx].get_c()[j] << endl;
        for (int k = 0; k < ChainNodeSPH[u].size(); k++){
            forchains.push_back(ChainNodeSPH[u][k]);
        }
    }
    for (auto it = forcycles.begin(); it != forcycles.end(); it++){
        AllccVars+= cyvar[*it];
    }
    for (auto it = forchains.begin(); it != forchains.end(); it++){
        AllccVars+= chvar[*it];
    }
    string name = "Active_CY_." + to_string(Cycles3rdTo2nd[idx] + 1);
    const char* cName = name.c_str();
    //cout << cyvar[Cycles3rdTo2nd[idx]].getName()  << ">=" <<  1 - ExprArcsVtx - AllccVars << endl;
    ActiveCCSubP_CY.add(IloRange(env, 1, cyvar[Cycles3rdTo2nd[idx]] + ExprArcsVtx + AllccVars, IloInfinity, cName));
    GrandSubProb.add(ActiveCCSubP_CY[ActiveCCSubP_CY.getSize() - 1]);
    ExprArcsVtx.end();
    AllccVars.end();
    
}
void Problem::GetNewIloRangeCH3rd(int idx, vector<Chain>&Chains2ndStage){
    IloExpr Expr2ArcsVtx (env, 0);
    IloExpr All2ccVars (env, 0);
    vector<int>forcycles;
    vector<int>forchains;
    for (int j = 0; j < Chains2ndStage[idx].Vnodes.size(); j++){
        int u = Chains2ndStage[idx].Vnodes[j].vertex;
        Expr2ArcsVtx+= vertex[u];
        if (j <= Chains2ndStage[idx].Vnodes.size() - 2){
            int v = Chains2ndStage[idx].Vnodes[j + 1].vertex;
            int s = mapArcs[make_pair(u,v)];
            Expr2ArcsVtx+= arc[u][s];
        }
        for (int k = 0; k < CycleNodeSPH[u].size(); k++){
            forcycles.push_back(CycleNodeSPH[u][k]);
        }
        for (int k = 0; k < ChainNodeSPH[u].size(); k++){
            if (ChainNodeSPH[u][k] != Chains3rdTo2nd[idx]) forchains.push_back(ChainNodeSPH[u][k]);
        }
    }
    for (auto it = forcycles.begin(); it != forcycles.end(); it++){
        All2ccVars+= cyvar[*it];
    }
    for (auto it = forchains.begin(); it != forchains.end(); it++){
        All2ccVars+= chvar[*it];
    }
    string name = "Active_CH_." + to_string(Chains3rdTo2nd[idx] + 1);
    const char* cName = name.c_str();
    //cout << chvar[Chains3rdTo2nd[idx]].getName()  << ">=" <<  1 - Expr2ArcsVtx - All2ccVars << endl;
    ActiveCCSubP_CH.add(IloRange(env, 1, chvar[Chains3rdTo2nd[idx]] + Expr2ArcsVtx + All2ccVars, IloInfinity, cName));
    GrandSubProb.add(ActiveCCSubP_CH[ActiveCCSubP_CH.getSize() - 1]);
    Expr2ArcsVtx.end();
    All2ccVars.end();
}
void Problem::GetNoGoodCut(map<pair<int,int>, bool>& FailedArcs, map<int, bool>& FailedVertices){
    IloExpr expr(env);
    
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        map<int, bool>::iterator it = FailedVertices.find(i);
        if (it == FailedVertices.end()){
            expr+= vertex[i];
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                map<pair<int,int>, bool>::iterator it2 = FailedArcs.find(make_pair(i, AdjacencyList[i][j] - 1));
                if (it2 == FailedArcs.end()) expr+= arc[i][j];
            }
        }
    }
    cout << expr << endl;
    string name = "NewBound." + to_string(ConsBeta.getSize());
    const char* cName = name.c_str();
    ConsBeta.add(IloRange(env, 0, Beta + (Nodes + NumArcs)*expr - TPMIP_Obj, IloInfinity, cName));
}
