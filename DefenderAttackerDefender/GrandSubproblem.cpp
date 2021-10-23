//
//  GrandSubproblem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-26.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandSubproblem.hpp"
void Problem::GrandSubProbMaster(){
    // Create model
    GrandSubProb = IloModel(env);
    cplexGrandSubP = IloCplex(GrandSubProb);
    cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, 1800);
    cplexGrandSubP.setParam(IloCplex::Param::Threads, 1);
    //Fill RepairedListCCs
    //FillRepSol_AND_CycleNodeGSP(CycleNodeGSP);
    FillRobustSolTHP();
    
    //Create cycle variables
    r = IloNumVarArray(env, GrandSubSolSet.size(), 0, 1, ILOINT);
    for (int i = 0; i < GrandSubSolSet.size(); i++){
       SetName1Index(r[i], "r", i);
       //cout << r[i].getName() << endl;
    }
    //Create arc variables
    arc = IloNumVarArray2(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        arc[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, 1, ILOINT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            //SetName2Index(arc[i][j], "arc",i, AdjacencyList[i][j] - 1);
        }
    }
    //Create vertex variables
    vertex = IloNumVarArray(env, Pairs, 0, 1, ILOINT);
    for (int i = 0; i < Pairs; i++){
        SetName1Index( vertex[i], "vertex",i);
    }
    //Create bounding variable
    Beta = IloNumVar(env);
    SetName1Index(Beta, "Beta", 0);
    
    //Add Bounding Constraint
    //BoundObjGrandSubP = IloRangeArray(env);
    string name = "Obj_GrandSubP";
    const char* cName = name.c_str();
    IloExpr obj (env,0);
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        obj += GrandSubSolSet[i].get_w()*r[i];
    }
    ObjGrandSubP = IloObjective(env, obj, IloObjective::Minimize, cName);
    GrandSubProb.add(ObjGrandSubP);
    
    MakeOneFailGrandSubP = IloRangeArray(env);
    for (int ite = 1; ite <= RepSolCounter; ite++){
        IloExpr ArcVerSum (env, 0);
        //ArcVerSum = GenerateMakeOneFailConstraint(ite);
        //cout << ArcVerSum << endl;
        name = "MakeOneFail." + to_string(ite);
        cName = name.c_str();
        MakeOneFailGrandSubP.add(IloRange(env, 1, ArcVerSum, IloInfinity, cName));
    }
    GrandSubProb.add(MakeOneFailGrandSubP);
    
    //Select a cycle if has not failed and no other was selected
    ActiveCCSubP_LB = IloRangeArray(env, GrandSubSolSet.size());
    map<int,vector<int>>::iterator it;
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        IloExpr ExprArcsVtx (env, 0);
        //ExprArcsVtx = GenerateAllArcsVertices(i);
        IloExpr AllccVars (env, 0);
        map<int,int>WhichCCs;
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            //ccVarsOneVtx (i, GrandSubSolSet[i].get_cc()[j], WhichCCs);
        }
        //Generate Sum with CCs containing the vertices in the i-th CC
        //AllccVars = GenerateAllCCsVars(WhichCCs);
        //cout << AllccVars << endl;
        name = "ActiveCC_LB." + to_string(i);
        cName = name.c_str();
        ActiveCCSubP_LB[i] = IloRange(env, 1, r[i] + ExprArcsVtx + AllccVars, IloInfinity, cName);
    }
    GrandSubProb.add(ActiveCCSubP_LB);
    
    //Guarantee one cycle per node
    float n;
    ActiveCCSubP_UB = IloRangeArray(env, GrandSubSolSet.size());
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        if (MaxArcFailures + MaxVertexFailures <= GrandSubSolSet[i].get_cc().size()*2){
            n = MaxArcFailures + MaxVertexFailures;
        }
        else{
            n = GrandSubSolSet[i].get_cc().size()*2;
        }
        IloExpr ExprArcsVtx (env, 0);
        IloExpr Total (env, 0);
        //ExprArcsVtx += GenerateAllArcsVertices(i);
        Total += ExprArcsVtx/n;
        name = "ActiveCC_UB." + to_string(i);
        cName = name.c_str();
        ActiveCCSubP_UB[i] = IloRange(env, -IloInfinity, r[i] + Total, 1, cName);
    }
    GrandSubProb.add(ActiveCCSubP_UB);
    
    //One cycle per vertex
    TheOneCC = IloRangeArray(env, Pairs);
    for (int i = 0; i < Pairs; i++){
        it = CycleNodeGSP.find(i);
        name = "TheOneCC." + to_string(i);
        cName = name.c_str();
        IloExpr Expr (env, 0);
        if (it != CycleNodeGSP.end()){
            for (int j = 0; j < CycleNodeGSP[i].size(); j++){
                Expr += r[CycleNodeGSP[i][j]];
            }
            TheOneCC[i] = IloRange(env, -IloInfinity, Expr, 1, cName);
        }
        else{
            TheOneCC[i] = IloRange(env, -IloInfinity, Expr, 1, cName);
        }
    }
    GrandSubProb.add(TheOneCC);
    
    IloExpr sumVertices (env, 0);
    for (int i = 0; i < Pairs; i++){
        sumVertices+=  vertex[i];
    }
    GrandSubProb.add(sumVertices <= MaxVertexFailures);
    
    IloExpr sumArcs (env, 0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumArcs+= arc[i][j];
        }
    }
    GrandSubProb.add(sumArcs <= MaxArcFailures);
    
    
    for (int i = 0; i < Pairs; i++){
        IloExpr sumDisjunctive (env, 0);
        //Outside
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumDisjunctive+= arc[i][j];
        }
        //Inside
        for (int j = 0; j < PredList[i].size(); j++){
            //int pos = FindPosArray(AdjacencyList[PredList[i][j]], i + 1);
            //sumDisjunctive+= arc[PredList[i][j]][pos];
        }
        name = "Disjunctive." + to_string(i);
        cName = name.c_str();
        GrandSubProb.add(IloRange(env, -IloInfinity, sumDisjunctive - (1 - vertex[i])*(AdjacencyList[i].getSize() + PredList[i].size()), 0, cName));
    }
    
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    IloNum GrandSubObj = 0;
    while(true){
        cplexGrandSubP.solve();
        if (cplexGrandP.getStatus() == 'Infeasible'){
            break;
        }
        else{
            GrandSubObj = cplexGrandSubP.getObjValue();
            env.out() << "Objective: " << GrandSubObj << endl;
            
            //Retrieve solution
            IloNumArray r_sol(env, GrandSubSolSet.size());
            cplexGrandSubP.getValues(r_sol,r);
            
            IloNumArray vertex_sol(env, Pairs);
            cplexGrandSubP.getValues(vertex_sol,vertex);
            
            IloNumArray2 arc_sol(env, AdjacencyList.getSize());
            for (int f = 0; f < arc_sol.getSize(); f++){
                arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
                cplexGrandSubP.getValues(arc_sol[f],arc[f]);
            }
            //Print curret solution
            cplexGrandSubP.exportModel("GrandSubP.lp");
            //Save new solution
            UpdateSNPSol(r_sol, GrandSubObj);
            PrintSolSNP(vertex_sol, arc_sol);

            
            RepSolCounter++;
            //BackRecoursePolicy(vertex_sol, arc_sol);
            //AllPolicy(vertex_sol, arc_sol);
            //RepairedListCCs = AmongPolicy(arc_sol); // It updates NewListCCs
            //Update arrays in GranSubSolSet and CycleNodeGSP
            //CheckNewIncludedVerticesGSP(GrandSubSolSet, RepairedListCCs);
            
            //Add new cols and rows
            //AddRow_MakeOneFail();
            //AddRowsCols_ActiveCCSubP_LB();
            cplexGrandSubP.exportModel("GrandSubP.lp");
        }
    }
    cout << "Hola" << endl;
    
    
}
void Problem::FillRobustSolTHP(){
    for(int i = 0; i < GrandSubSolSet.size(); i++){
        RobustSolTHP.push_back(Cycles(GrandSubSolSet[i].get_cc(), GrandSubSolSet[i].get_w()));
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
            RobustSolTHP.push_back(Cycles(GrandSubSolSet[i].get_cc(), GrandSubSolSet[i].get_w()));
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
vector<Cycles> Problem::BackRecoursePolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    
    
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
                NewListCCs.back().set_Weight(int(NewList[k].get_c().size()));
                for (int l = 0; l < NewList[k].get_c().size(); l++){
                    CycleNodeSPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
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
    
    IloNumArray2 AdjaList (env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env);
    }
    
    
    //Fill in Adjacency List
    for (int i = 0; i < ListVertices.size(); i++){
        //cout << endl << LookFor[i] << "\t" << ":";
        for (int j = 0; j < AdjacencyList[ListVertices[i]].getSize(); j++){
            AdjaList[ListVertices[i]].add(AdjacencyList[ListVertices[i]][j]);
        }
    }
    
    //Find among-cycles
    vector<Cycles> NewList;
    for (int i = 0; i < ListVertices.size(); i++){
        int origin = ListVertices[i];
        NewList = SubCycleFinder(env, AdjaList, origin);
        for (int k = 0; k < NewList.size(); k++){
            NewListCCs.push_back(NewList[k]);
            NewListCCs.back().set_Weight(int(NewList[k].get_c().size()));
            for (int l = 0; l < NewList[k].get_c().size(); l++){
                CycleNodeSPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
            }
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    
    return NewListCCs;
}
vector<Cycles> Problem::AllPolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;

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
                CycleNodeSPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
            }
            NewListCCs.back().set_Weight(counter);
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    return NewListCCs;
}
void Problem::InitializeVertexinSolChain(vector<int>&ListVertices,vector<vChain>& VertexinSolChain){
    
    for (int j = 0; j < ListVertices.size(); j++){
        if (VertexinSolChain[ListVertices[j]].veci.size() > 0){
            VertexinSolChain[ListVertices[j]].veci = vector<int>();
        }
        for (int l = 0; l < AdjacencyList[ListVertices[j]].getSize(); l++){
            int neighbour = AdjacencyList[ListVertices[j]][l] - 1;
            int pos = FindPosVector(ListVertices, neighbour);
            if (pos != -1) {
                VertexinSolChain[ListVertices[j]].veci.push_back(neighbour);
            }
        }
        VertexinSolChain[ListVertices[j]].it = VertexinSolChain[ListVertices[j]].veci.begin();
        VertexinSolChain[ListVertices[j]].itEnd = VertexinSolChain[ListVertices[j]].veci.end();
    }
}
vector<vector<int>> Problem::FindChains(IloNumArray2 x_sol, vector<vChain>& VertexinSolChain, vector<int>& vinFirstStageSol){
    vector<vector<int>>Bestchains;
    int u,v;
    double w = 0, price = 0;
    vector<vChain>::iterator itVinChain = VertexinSolChain.begin();
    vector<vChain>::iterator itVinChainAux;
    vector<Chain>PPChains;
    
    while (itVinChain->vertex > Pairs - 1){//By construction altruistic donors come first
        if (itVinChain->veci.size() > 0){
            PPChains.push_back(Chain(*itVinChain));
            PPChains.back().AccumWeight++;
            while(PPChains.back().Vnodes.size() >  0){//next neighbor
                //Increase iterator
                if (PPChains.back().Vnodes.back().FirstPass == true){
                    PPChains.back().Vnodes.back().FirstPass = false;
                }else{PPChains.back().Vnodes.back().it++;}
                u = PPChains.back().Vnodes.back().vertex;
                v = *(PPChains.back().Vnodes.back().it);
                w = Weights[make_pair(u, v)];
                //Find pointer to v in VertexinSolChain
                itVinChainAux = VertexinSolChain.begin();
                for (itVinChainAux; itVinChainAux != VertexinSolChain.end(); itVinChainAux++)
                    if (itVinChainAux->vertex == v){
                        break;
                    }
                //If chain is not too long or we have reached the end of some vertex's neighbors
                if (PPChains.back().Vnodes.size() + 1 <= ChainLength  +  1 || PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().veci.end()){
                    //If vertex not already in chain
                    if (v2AlreadyinChain(PPChains.back().Vnodes, v) == false){
                    //If vertex v is in VertexinSolChain add vertex to current chain, augment AccumWeight, store current chain
                        //Add vertex to current chain
                        PPChains.back().Vnodes.push_back(*itVinChainAux);
                        if (v2inFirstStageSol(vinFirstStageSol, v) == true){// Update weight only if v is in the 1st-stage solution
                            PPChains.back().AccumWeight++;
                        }
                    }
                }
                else{
                    PPChains.back().Vnodes.pop_back();
                    //Find new neighbor
                    FindNewNeighbor(PPChains);
                }
            }
            PPChains.erase(PPChains.end() - 1);
        }
        itVinChain++;
    }
    //Transfer chains to Bestchains
    for (int i = 0; i < PPChains.size();i++){
        Bestchains.push_back(vector<int>());
        for (int j = 0; j < PPChains[i].Vnodes.size(); j++){
            Bestchains.back().push_back(PPChains[i].Vnodes[j].vertex);
        }
    }
    
    if (Bestchains.size() > 0){
        if (Bestchains[0].size() == 1){
            cout << "S.O.S" << endl;
        }
    }
    return Bestchains;
    
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
