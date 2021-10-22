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
vector<Cycles> Problem::BackRecoursePolicy(vector<vector<int>>& FistStageSol){
    vector<Cycles> NewListCCs;
    
    for (int i = 0; i < FistStageSol.size(); i++){
        IloNumArray2 AdjaList (env, AdjacencyList.getSize());

        for (int i = 0; i < AdjacencyList.getSize(); i++){
            AdjaList[i] = IloNumArray(env);
        }
        for (int j = 0; j < FistStageSol[i].size(); j++){
            for (int l = 0; l < AdjacencyList[FistStageSol[i][j]].getSize(); l++){
                int neighbour = AdjacencyList[FistStageSol[i][j]][l] - 1;
                int pos = FindPosVector(FistStageSol[i], neighbour);
                if (pos != -1) {
                    AdjaList[FistStageSol[i][j]].add(neighbour + 1);
                }
            }
        }
//        for (int i = 0; i < AdjaList.getSize(); i++){
//            cout << endl << i << ": ";
//            for (int j = 0; j < AdjaList[i].getSize(); j++){
//                cout << AdjaList[i][j] << "\t";
//            }
//        }
        
        //Find inner cycles
        vector<Cycles> NewList;
        for (int j = 0; j < FistStageSol[i].size(); j++){
            int origin = FistStageSol[i][j];
            NewList = SubCycleFinder(env, AdjaList, origin);
            for (int k = 0; k < NewList.size(); k++){
                NewListCCs.push_back(NewList[k]);
                NewListCCs.back().set_Many(int(NewList[k].get_c().size()));
                for (int l = 0; l < NewList[k].get_c().size(); l++){
                    CycleNodeTPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
                }
            }
            //Remove origin from AdjaList
            AdjaList[origin].clear();
        }
        AdjaList.clear();
    }
    
//    for (int k = 0; k < NewListCCs.size(); k++){
//        cout << endl;
//        for (int l = 0; l < NewListCCs[k].get_c().size(); l++){
//            cout << NewListCCs[k].get_c()[l] << "\t";
//        }
//    }
    return NewListCCs;
}
vector<Cycles> Problem::AmongPolicy(vector<vector<int>>& FistStageSol){
    vector<Cycles> NewListCCs;
    
    IloNumArray2 AdjaList (env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env);
    }
    
    vector<int>LookFor;
    
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            LookFor.push_back(GrandSubSolSet[i].get_cc()[j]);
        }
    }
    
    //Fill in Adjacency List
    for (int i = 0; i < LookFor.size(); i++){
        //cout << endl << LookFor[i] << "\t" << ":";
        for (int j = 0; j < AdjacencyList[LookFor[i]].getSize(); j++){
            int pos = FindPosVector(LookFor, AdjacencyList[LookFor[i]][j] - 1);
            if (pos != -1){
                AdjaList[LookFor[i]].add(AdjacencyList[LookFor[i]][j]);
            }
        }
    }
    
    //Find among-cycles
    vector<Cycles> NewList;
    for (int i = 0; i < LookFor.size(); i++){
        int origin = LookFor[i];
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
vector<Cycles> Problem::AllPolicy(vector<vector<int>>& FistStageSol){
    vector<Cycles> NewListCCs;
    vector<int>ListVertices;
    
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            ListVertices.push_back(GrandSubSolSet[i].get_cc()[j]);
        }
    }
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
                CycleNodeTPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
            }
            NewListCCs.back().set_Many(counter);
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    return NewListCCs;
}
