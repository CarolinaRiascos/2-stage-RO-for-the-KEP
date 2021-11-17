//
//  GrandSubproblem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-26.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandSubproblem.hpp"
void Problem::GrandSubProbMaster(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<IndexGrandSubSol>&SolFirstStage){
    // Create model
    GrandSubProb = IloModel(env);
    cplexGrandSubP = IloCplex(GrandSubProb);
    cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, 1800);
    cplexGrandSubP.setParam(IloCplex::Param::Threads, 1);
    
    SampleCols2ndStage(Chains2ndStage, Cycles2ndStage);
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
    arc = IloNumVarArray2(env,AdjacencyList.getSize());
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
    //Create selected-vertex variables
    vector<int>ListSelVertices = GetSelVertices(SolFirstStage);
    selvertex = IloNumVarArray(env, ListSelVertices.size(), 0, 1, ILOINT);
    for (int i = 0; i < selvertex.getSize(); i++){
        SetName1Index(selvertex[i], "selv", ListSelVertices[i] + 1);
    }
    
    //Create bounding variable
    Beta = IloNumVar(env);
    SetName1Index(Beta, "Beta", 0);
    //Create excess variable in objective
    Excess = IloNumVar(env);
    IloExpr excess (env, 0);
    SetName1Index(Excess, "Excess", 0);
    for (int i = 0; i < ListSelVertices.size(); i++){
        excess += 1 - selvertex[i];
    }
    GrandSubProb.add(IloRange(env, 0, Excess - excess, 0));
    SelVert2ndPH = IloRangeArray(env, 0);
    //Add vertex-selection constraint
    for (int i = 0; i < ListSelVertices.size(); i++){
        IloExpr exSel(env, 0);
        for (int j = 0; j < CycleNodeSPH[ListSelVertices[i]].size(); j++){
            exSel+= cyvar[CycleNodeSPH[ListSelVertices[i]][j]];
        }
        for (int j = 0; j < ChainNodeSPH[ListSelVertices[i]].size(); j++){
            exSel+= chvar[ChainNodeSPH[ListSelVertices[i]][j]];
        }
        SelVert2ndPH.add(IloRange(env, -IloInfinity, selvertex[i] - exSel, 0));
    }
    GrandSubProb.add(SelVert2ndPH);
    //Add Objective function and Bounding Constraint
    string name = "Obj_GrandSubP";
    
    const char* cName = name.c_str();
    IloExpr obj (env,0);
    for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
        obj += Cycles2ndStage[Cycles2ndTo3rd[i]].get_Many()*cyvar[i];
    }
    for (int i = 0; i < Chains2ndTo3rd.size(); i++){
        obj += Chains2ndStage[Chains2ndTo3rd[i]].AccumWeight*chvar[i];
    }
    ObjGrandSubP = IloObjective(env, Beta + Excess, IloObjective::Minimize, cName);
    GrandSubProb.add(ObjGrandSubP);
    ConsBeta.add(IloRange(env, 0, Beta - obj , IloInfinity));
    GrandSubProb.add(ConsBeta);
    //obj.end();
    
    //Select a cycle/chain if it has not failed and no other was selected
    ActiveCCSubP_CY = IloRangeArray(env);
    ActiveCCSubP_CH = IloRangeArray(env);
    for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
        IloExpr ExprArcsVtx (env, 0);
        IloExpr AllccVars (env, 0);
        map<int,bool>forcycles;
        map<int,bool>forchains;
        int idx = Cycles2ndTo3rd[i];
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
                if (CycleNodeSPH[u][k] != i) forcycles[CycleNodeSPH[u][k]];
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
        name = "Active_cyvar2nd." + to_string(i);
        cName = name.c_str();
        //cout << cyvar[i].getName()  << ">=" <<  1 - ExprArcsVtx - AllccVars << endl;
        ActiveCCSubP_CY.add(IloRange(env, 1, cyvar[i] + ExprArcsVtx + AllccVars, IloInfinity, cName));
        GrandSubProb.add(ActiveCCSubP_CY[ActiveCCSubP_CY.getSize() - 1]);
    }
    //Chain variables
    for (int i = 0; i < Chains2ndTo3rd.size(); i++){
        IloExpr Expr2ArcsVtx (env, 0);
        IloExpr All2ccVars (env, 0);
        map<int,bool>forcycles;
        map<int,bool>forchains;
        int idx = Chains2ndTo3rd[i];
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
                if (ChainNodeSPH[u][k] != i) forchains[ChainNodeSPH[u][k]];
            }
        }
        for (auto it = forcycles.begin(); it != forcycles.end(); it++){
            All2ccVars+= cyvar[it->first];
        }
        for (auto it = forchains.begin(); it != forchains.end(); it++){
            All2ccVars+= chvar[it->first];
        }
        name = "Active_chvar2nd." + to_string(i);
        cName = name.c_str();
        //cout << cyvar[i].getName()  << ">=" <<  1 - Expr2ArcsVtx - All2ccVars << endl;
        ActiveCCSubP_CH.add(IloRange(env, 1, chvar[i] + Expr2ArcsVtx + All2ccVars, IloInfinity, cName));
        GrandSubProb.add(ActiveCCSubP_CH[ActiveCCSubP_CH.getSize() - 1]);
    }
    
    IloExpr sumVertices (env, 0);
    for (int i = 0; i < Nodes; i++){
        sumVertices+=  vertex[i];
    }
    GrandSubProb.add(sumVertices <= 0);
    
    IloExpr sumArcs (env, 0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumArcs+= arc[i][j];
        }
    }
    GrandSubProb.add(sumArcs <= 0);
    
    cplexGrandSubP.exportModel("GrandSubP.lp");
    cplexGrandSubP.solve();
    if (cplexGrandSubP.getStatus() == 'Infeasible'){
        cout << "S.O.S. This should not happen." << endl;
    }
    else{
        SPMIP_Obj = cplexGrandSubP.getValue(Beta);
        env.out() << "Objective: " << SPMIP_Obj << endl;//cplexGrandSubP.getValue(Excess)
        
        //Retrieve solution
        IloNumArray cyvar_sol(env, Cycles2ndTo3rd.size());
        IloNumArray chvar_sol(env, Chains2ndTo3rd.size());
        cplexGrandSubP.getValues(cyvar_sol,cyvar);
        cplexGrandSubP.getValues(chvar_sol,chvar);
        
        for (int i = 0; i < cyvar_sol.getSize(); i++){
            if (cyvar_sol[i] > 0.9) cout << cyvar[i].getName() << endl;
        }
        for (int i = 0; i < chvar_sol.getSize(); i++){
            if (chvar_sol[i] > 0.9) cout << chvar[i].getName() << endl;
        }
        cout << "Beta: " << cplexGrandSubP.getValue(Beta) << endl;
        
        IloNumArray vertex_sol(env, Pairs);
        cplexGrandSubP.getValues(vertex_sol,vertex);
        
        IloNumArray2 arc_sol(env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
            cplexGrandSubP.getValues(arc_sol[f],arc[f]);
        }

        RepSolCounter++;
    }
    cout << "Hola" << endl;
}
void Problem::FillRobustSolTHP(){
    for(int i = 0; i < GrandProbSol.size(); i++){
        //RobustSolTHP.push_back(Cycles(GrandProbSol[i].get_cc(), GrandProbSol[i].get_w()));
    }
}
void Problem::PrintSolSNP(IloNumArray vertex_sol, IloNumArray2 arc_sol){
    cout << "Iteration: " << RepSolCounter << endl;
    for (int i = 0; i < RobustSolTHP.size(); i++){
        cout << "c." << i << " :  ";
        for (int j = 0; j < RobustSolTHP[i].get_c().size(); j++){
            cout << RobustSolTHP[i].get_c()[j] << "\t";
        }
        cout << endl;
    }
    cout << "Vertices failed: ";
    for (int i = 0; i < vertex_sol.getSize(); i++){
        if (vertex_sol[i] > 0.9){
            cout << i << endl;
        }
    }
    cout << "Arcs failed: ";
    for (int i = 0; i < arc_sol.getSize(); i++){
        for (int j = 0; j < arc_sol[i].getSize(); j++){
            if (arc_sol[i][j] > 0.9){
                cout << i << "->" << AdjacencyList[i][j] - 1 << "\t";
            }
        }
    }
    
}
void Problem::UpdateSNPSol(IloNumArray& r_sol, IloNum GrandSubObj){
    RobustObjTPH = GrandSubObj;
    RobustSolTHP.clear();
    for (int i = 0; i < r_sol.getSize(); i++){
        if (r_sol[i] > 0.9){
            //RobustSolTHP.push_back(Cycles(GrandProbSol[i].get_cc(), GrandProbSol[i].get_w()));
        }
    }
}
int FindPosVector(vector<int> array, int value){
    int pos = -1;
    for (int j = 0; j < array.size(); j++){
        if (array[j] == value) return j;
    }
    return pos;
}
vector<Cycles> Problem::Get2ndStageCycles (vector<IndexGrandSubSol>& GrandProbSol, string policy){
    //Create List of all vertices in 1st-stage sol: vinFirstStageSol
    vector<int> vinFirstStageSol;
    if (policy == "among" || policy == "full"){
        for (int i = 0; i < GrandProbSol.size(); i++){
            for (int j = 0; j < GrandProbSol[i].get_cc().size(); j++){
                vinFirstStageSol.push_back(GrandProbSol[i].get_cc()[j]);
            }
        }
    }
    
    //Create List of all candidate vertices: ListVertices
    vector<int>ListVertices;
    if (policy == "among"){
        ListVertices = vinFirstStageSol;
    }
    else if (policy == "full"){
        for (int i = 0; i < Pairs; i++){
            ListVertices.push_back(i);
        }
    }
    
    vector<Cycles>RecoCycles;
    if (policy == "among"){
        //Call Find Cycles
        RecoCycles = AmongPolicy(vinFirstStageSol);
    }
    else if (policy == "full"){
        //Call Find Cycles
        RecoCycles = AllPolicy(vinFirstStageSol);
    }
    else{//BackArcs recourse
        for (int i = 0; i < GrandProbSol.size(); i++){
            //Call Find Cycles
            vector<Cycles>auxCycles;
            vector<int>aux = GrandProbSol[i].get_cc();
            auxCycles = BackRecoursePolicy(aux);
            for (int j = 0; j < auxCycles.size(); j++){
                RecoCycles.push_back(auxCycles[j]);
            }
        }
    }
    //Fill in CycleNodeTPH
    for (int i = 0; i < RecoCycles.size(); i++){
        for (int j = 0; j < RecoCycles[i].get_c().size(); j++){
            CycleNodeTPH[RecoCycles[i].get_c()[j]].push_back(i);
            if (j <= RecoCycles[i].get_c().size() - 2){
                ArcsinCyclesTHP[make_pair(RecoCycles[i].get_c()[j], RecoCycles[i].get_c()[j + 1])].push_back(i);
            }
            else{
                ArcsinCyclesTHP[make_pair(RecoCycles[i].get_c()[j], RecoCycles[i].get_c()[0])].push_back(i);
            }
        }
    }
    
    return RecoCycles;
}
vector<Cycles> Problem::BackRecoursePolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();
    
        IloNumArray2 AdjaList (env, AdjacencyList.getSize());

        for (int j = 0; j < AdjacencyList.getSize(); j++){
            AdjaList[j] = IloNumArray(env);
        }
        for (int j = 0; j < ListVertices.size(); j++){
            for (int l = 0; l < AdjacencyList[ListVertices[j]].getSize(); l++){
                int neighbour = AdjacencyList[ListVertices[j]][l] - 1;
                int pos = FindPosVector(ListVertices, neighbour);
                if (pos != -1) {
                    AdjaList[ListVertices[j]].add(neighbour + 1);
                }
            }
        }
        
        //Find inner cycles
        vector<Cycles> NewList;
        for (int j = 0; j < ListVertices.size(); j++){
            int origin = ListVertices[j];
            NewList = SubCycleFinder(env, AdjaList, origin);
            for (int k = 0; k < NewList.size(); k++){
                NewListCCs.push_back(NewList[k]);
                NewListCCs.back().set_Many(int(NewList[k].get_c().size()));
                for (int l = 0; l < NewList[k].get_c().size(); l++){
                    CycleNodeTPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
                    if (l <= NewList[k].get_c().size() - 2){
                        ArcsinCyclesTHP[make_pair(NewList[k].get_c()[l], NewList[k].get_c()[l + 1])].push_back(int(NewListCCs.size() - 1));
                    }
                    else{
                        ArcsinCyclesTHP[make_pair(NewList[k].get_c()[l], NewList[k].get_c()[0])].push_back(int(NewListCCs.size() - 1));
                    }
                }
            }
            //Remove origin from AdjaList
            AdjaList[origin].clear();
        }
        AdjaList.clear();
    
    return NewListCCs;
}
vector<Cycles> Problem::AmongPolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();
    
    IloNumArray2 AdjaList (env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env);
    }
    
    
    //Fill in Adjacency List
    for (int i = 0; i < ListVertices.size(); i++){
        //cout << endl << LookFor[i] << "\t" << ":";
        for (int j = 0; j < AdjacencyList[ListVertices[i]].getSize(); j++){
            int pos = FindPosVector(ListVertices, AdjacencyList[ListVertices[i]][j] - 1);
            if (pos != -1){
                AdjaList[ListVertices[i]].add(AdjacencyList[ListVertices[i]][j]);
            }
        }
    }
    
    //Find among-cycles
    vector<Cycles> NewList;
    for (int i = 0; i < ListVertices.size(); i++){
        int origin = ListVertices[i];
        NewList = SubCycleFinder(env, AdjaList, origin);
        for (int k = 0; k < NewList.size(); k++){
            NewListCCs.push_back(NewList[k]);
            NewListCCs.back().set_Many(int(NewList[k].get_c().size()));
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    
    return NewListCCs;
}
vector<Cycles> Problem::AllPolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();

    //Make a copy of Adjacency List
    IloNumArray2 AdjaList(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env, AdjacencyList[i].getSize());
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            AdjaList[i][j] = AdjacencyList[i][j];
        }
    }
    
    //Find all cycles
    vector<Cycles> NewList;
    for (int i = 0; i < ListVertices.size(); i++){
        int origin = ListVertices[i];
        NewList = SubCycleFinder(env, AdjaList, origin);
        for (int k = 0; k < NewList.size(); k++){
            NewListCCs.push_back(NewList[k]);
            int counter = 0;
            for (int l = 0; l < NewList[k].get_c().size(); l++){
                int pos = FindPosVector(ListVertices, NewList[k].get_c()[l]);
                if (pos != -1) counter++;
            }
            NewListCCs.back().set_Many(counter);
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    return NewListCCs;
}
vector<Chain> Problem::Get2ndStageChains (vector<IndexGrandSubSol>& GrandProbSol, string policy){
    
    //Create List of all vertices in 1st-stage sol: vinFirstStageSol
    vector<int> vinFirstStageSol;
    if (policy == "among" || policy == "full"){
        for (int i = 0; i < GrandProbSol.size(); i++){
            for (int j = 0; j < GrandProbSol[i].get_cc().size(); j++){
                vinFirstStageSol.push_back(GrandProbSol[i].get_cc()[j]);
            }
        }
    }
    
    //Create List of all candidate vertices: ListVertices
    vector<int>ListVertices;
    if (policy == "among"){
        ListVertices = vinFirstStageSol;
    }
    else if (policy == "full"){
        for (int i = 0; i < Nodes; i++){
            ListVertices.push_back(i);
        }
    }
    
    VertexinSolChain = vector<vChain>();
    vector<Chain>RecoChains;
    if (policy == "among" || policy == "full"){
        //Call InitializeVertexinSolChain
        InitializeVertexinSolChain(ListVertices, VertexinSolChain);
        //Call Find Chains
        RecoChains = FindChains(VertexinSolChain, vinFirstStageSol);
    }else{//BackArcs recourse
        for (int i = 0; i < GrandProbSol.size(); i++){
            //Call InitializeVertexinSolChain
            vector<int> aux = GrandProbSol[i].get_cc();
            InitializeVertexinSolChain(aux, VertexinSolChain);
            //Call Find Chains
            vector<Chain>auxChains;
            auxChains = FindChains(VertexinSolChain, aux);
            for (int j = 0; j < auxChains.size(); j++){
                RecoChains.push_back(auxChains[j]);
            }
        }
    }
    
    //Fill in ChainNodeTPH
    for (int i = 0; i < RecoChains.size(); i++){
        for (int j = 0; j < RecoChains[i].Vnodes.size(); j++){
            ChainNodeTPH[RecoChains[i].Vnodes[j].vertex].push_back(i);
            if (j <= RecoChains[i].Vnodes.size() - 2){
                ArcsinChainsTHP[make_pair(RecoChains[i].Vnodes[j].vertex, RecoChains[i].Vnodes[j + 1].vertex)].push_back(i);
            }
            else{
                ArcsinChainsTHP[make_pair(RecoChains[i].Vnodes[j].vertex, RecoChains[i].Vnodes[0].vertex)].push_back(i);
            }
        }
    }
    
    return RecoChains;
}
void Problem::InitializeVertexinSolChain(vector<int>&ListVertices,vector<vChain>& VertexinSolChain){
    vector<int> null(1, -1);
    for (int j = 0; j < ListVertices.size(); j++){
        vector<int> veci;
        for (int l = 0; l < AdjacencyList[ListVertices[j]].getSize(); l++){
            int neighbour = AdjacencyList[ListVertices[j]][l] - 1;
            int pos = FindPosVector(ListVertices, neighbour);
            if (pos != -1) {
                veci.push_back(neighbour);
            }
        }
        
        if (veci.size() > 0){
            VertexinSolChain.push_back(vChain(ListVertices[j], veci));
            VertexinSolChain.back().it = VertexinSolChain[ListVertices[j]].veci.begin();
            VertexinSolChain.back().itEnd = VertexinSolChain[ListVertices[j]].veci.end();
        }
        else{
            VertexinSolChain.push_back(vChain(ListVertices[j], null));
            VertexinSolChain.back().it = VertexinSolChain[ListVertices[j]].veci.begin();
            VertexinSolChain.back().itEnd = VertexinSolChain[ListVertices[j]].veci.end();
        }
    }
}
vector<Chain> Problem::FindChains(vector<vChain>& VertexinSolChain, vector<int>& vinFirstStageSol){
    vector<vector<int>>Chains;
    vector<int>null(0);
    int u,v;
    double w = 0, price = 0; bool isin = false;
    vector<vChain>::iterator itVinChain = VertexinSolChain.begin();
    vector<vChain>::iterator itVinChainAux;
    vector<Chain>PPChains;
    ChainNodeTPH.clear();
    ArcsinChainsTHP.clear();
    bool justPopedBacked = false;
    
    for (itVinChain; itVinChain != VertexinSolChain.end(); itVinChain++){//By construction altruistic donors come first
        if (itVinChain->vertex > Pairs - 1 && v2inFirstStageSol(vinFirstStageSol, itVinChain->vertex) == true){
            PPChains.push_back(Chain(*itVinChain));
            while(PPChains.back().Vnodes.size() >  0){//next neighbor
                //Increase iterator
                if (PPChains.back().Vnodes.back().FirstPass == true){
                    PPChains.back().Vnodes.back().FirstPass = false;
                }else{PPChains.back().Vnodes.back().it++;}
                u = PPChains.back().Vnodes.back().vertex;
                v = *(PPChains.back().Vnodes.back().it);
                if (PPChains.back().Vnodes.back().veci.size() == 1 && v == -1){
                    PPChains.back().Vnodes.back().it = PPChains.back().Vnodes.back().itEnd;
                }
                //Find pointer to v in VertexinSolChain
                itVinChainAux = VertexinSolChain.begin();
                for (itVinChainAux; itVinChainAux != VertexinSolChain.end(); itVinChainAux++)
                    if (itVinChainAux->vertex == v){
                        break;
                    }
                //If chain is not too long or we have reached the end of some vertex's neighbors
                if (PPChains.back().Vnodes.size() + 1 <= ChainLength  +  1 && PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().itEnd){
                    //If vertex not already in chain
                    if (v2AlreadyinChain(PPChains.back().Vnodes, v) == false){
                    //If vertex v is in VertexinSolChain add vertex to current chain, augment AccumWeight, store current chain
                        isin = v2inFirstStageSol(vinFirstStageSol, v);
                        if (isin == false && PPChains.back().Vnodes.size() + 1 == ChainLength  +  1){
                            //Do nothing
                        }
                        else{
                            if (isin == true && PPChains.back().AccumWeight >= 1 && justPopedBacked == false){// Update weight only if v is in the 1st-stage solution
                                Chain aux = PPChains.back();
                                PPChains.push_back(aux);
                                PPChains.back().AccumWeight++;
                                for (int j = 0; j < PPChains.back().Vnodes.size(); j++){
                                    int v = PPChains.back().Vnodes[j].vertex;
                                    //Map arcs to cycles
                                    if (j <= PPChains.back().Vnodes.size() - 2){
                                        ArcsinChainsTHP[make_pair(v, PPChains.back().Vnodes[j + 1].vertex)].push_back(int(PPChains.size() - 1));
                                    }
                                    else{
                                        ArcsinChainsTHP[make_pair(v, itVinChainAux->vertex)].push_back(int(PPChains.size() - 1));
                                    }
                                }
                            }
                            else {
                                if (isin == true)   PPChains.back().AccumWeight++;
                            }
                            //Add vertex to current chain
                            PPChains.back().Vnodes.push_back(*itVinChainAux);
                            justPopedBacked = false;
                        }
                    }
                }
                else{
                    if (PPChains.back().AccumWeight > 0){
                        isin = v2inFirstStageSol(vinFirstStageSol, PPChains.back().Vnodes.back().vertex);
                        if (isin == false){
                            PPChains.back().Vnodes.pop_back();
                            Chain aux = PPChains.back();
                            if (justPopedBacked == false) PPChains.push_back(aux);
                        }
                        else{
                            Chain aux = PPChains.back();
                            if (justPopedBacked == false) PPChains.push_back(aux);
                            PPChains.back().Vnodes.pop_back();
                            PPChains.back().AccumWeight--;
                        }
                    }
                    else{
                        PPChains.back().Vnodes.pop_back();
                    }
                    justPopedBacked = true;
                //Find new neighbor
                //FindNewNeighbor(PPChains);
                }
            }
            PPChains.erase(PPChains.end() - 1);
        }
    }
    
    if (Chains.size() > 0){
        if (Chains[0].size() == 1){
            cout << "S.O.S" << endl;
        }
    }
    return PPChains;
}
bool v2inFirstStageSol(vector<int>sol, int v){
    for (int i = 0; i < sol.size(); i++){
        if (v == sol[i]) return true;
    }
    return false;
}
bool v2AlreadyinChain(vector<vChain> v1, int v2){
    for (int j = 0; j < v1.size(); j++){
        if (v1[j].vertex == v2){
            return true;
        }
    }
    return false;
}
void FindNewNeighbor(vector<Chain>& PPChains){
    while (PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().veci.end()){
        PPChains.back().Vnodes.back().it++;
        if (PPChains.back().Vnodes.back().it == PPChains.back().Vnodes.back().itEnd){
            PPChains.back().Vnodes.pop_back();
        }
        else{
            break;
        }
        if (PPChains.back().Vnodes.size() == 0) break;
    }
}
bool sortChains(Chain& c1, Chain& c2){
    return (c1.AccumWeight > c2.AccumWeight);
}
bool sortCycles(Cycles& c1, Cycles& c2){
    return (c1.get_Many() > c2.get_Many());
}
void Problem::SampleCols2ndStage(vector<Chain>&Chains, vector<Cycles>&Cycles){
    //Organize vector by AccumWeight and HowMany
    sort(Chains.begin(), Chains.end(), sortChains);
    sort(Cycles.begin(), Cycles.end(), sortCycles);
    
    CycleNodeSPH.size();
    ChainNodeSPH.size();
    
    //Select 100 cycles and 200 chains
    int counter = 0;
    for (int i = 0; i < Cycles.size(); i++){
        Cycles2ndTo3rd[counter] = i;
        Cycles3rdTo2nd[i] = counter;
        for (int j = 0; j < Cycles[i].get_c().size(); j++){
            CycleNodeSPH[Cycles[i].get_c()[j]].push_back(counter);
        }
        counter++;
        if (counter == 100) break;
    }
    counter = 0;
    for (int i = 0; i < Chains.size(); i++){
        Chains2ndTo3rd[counter] = i;
        Cycles3rdTo2nd[i] = counter;
        for (int j = 0; j < Chains[i].Vnodes.size(); j++){
            int u = Chains[i].Vnodes[j].vertex;
            ChainNodeSPH[u].push_back(counter);
        }
        counter++;
        if (counter == 100) break;
    }
}
vector<int>Problem::GetSelVertices(vector<IndexGrandSubSol>&SolFirstStage){
    vector<int>ListSelVertices;
    for (auto it = SolFirstStage.begin(); it != SolFirstStage.end(); it++){
        for (auto it2 = it->get_cc().begin(); it2 != it->get_cc().end(); it2++){
            ListSelVertices.push_back(*it2);
        }
    }
    
    return ListSelVertices;
}

