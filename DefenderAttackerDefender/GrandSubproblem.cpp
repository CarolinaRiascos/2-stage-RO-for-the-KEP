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
    FillRepSol_AND_CycleNodeGSP(CycleNodeGSP);
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
            SetName2Index(arc[i][j], "arc",i, AdjacencyList[i][j] - 1);
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
        ArcVerSum = GenerateMakeOneFailConstraint(ite);
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
        ExprArcsVtx = GenerateAllArcsVertices(i);
        IloExpr AllccVars (env, 0);
        map<int,int>WhichCCs;
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            ccVarsOneVtx (i, GrandSubSolSet[i].get_cc()[j], WhichCCs);
        }
        //Generate Sum with CCs containing the vertices in the i-th CC
        AllccVars = GenerateAllCCsVars(WhichCCs);
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
        ExprArcsVtx += GenerateAllArcsVertices(i);
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
            int pos = FindPosArray(AdjacencyList[PredList[i][j]], i + 1);
            sumDisjunctive+= arc[PredList[i][j]][pos];
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
            RepairedListCCs = AmongPolicy(vertex_sol, arc_sol); // It updates NewListCCs
            //Update arrays in GranSubSolSet and CycleNodeGSP
            CheckNewIncludedVerticesGSP(GrandSubSolSet, RepairedListCCs);
            
            //Add new cols and rows
            AddRow_MakeOneFail();
            AddRowsCols_ActiveCCSubP_LB();
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
void Problem::AddRowsCols_ActiveCCSubP_LB(){
    vector<int> PosNewCycles;
    string name;
    const char* cName;
    //Add New Constraints ActiveCCSubP_LB
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        if (GrandSubSolSet[i].get_itev().size() == 1 && GrandSubSolSet[i].get_ite(0) == RepSolCounter){
            PosNewCycles.push_back(i);
            IloExpr ExprArcsVtx (env, 0);
            ExprArcsVtx = GenerateAllArcsVertices(i);
            IloExpr AllccVars (env, 0);
            map<int,int>WhichCCs;
            for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
                ccVarsOneVtx (i, GrandSubSolSet[i].get_cc()[j], WhichCCs);
            }
            //Generate Sum with CCs containing the vertices in the i-th CC
            AllccVars = GenerateAllCCsVars(WhichCCs);
            name = "ActiveCC_LB." + to_string(i);
            cName = name.c_str();
            ActiveCCSubP_LB.add(IloRange(env, 1, ExprArcsVtx + AllccVars, IloInfinity, cName));
            GrandSubProb.add(ActiveCCSubP_LB[ActiveCCSubP_LB.getSize() - 1]);
            //Geberate  ActiveCCSubP_UB
            int n = 0;
            if (MaxArcFailures + MaxVertexFailures <= GrandSubSolSet[i].get_cc().size()*2){
                n = int(MaxArcFailures + MaxVertexFailures);
            }
            else{
                n = int(GrandSubSolSet[i].get_cc().size())*2;
            }
            IloExpr Total (env, 0);
            Total += ExprArcsVtx/n;
            //cout << Total << endl;
            name = "ActiveCC_UB." + to_string(i);
            cName = name.c_str();
            ActiveCCSubP_UB.add(IloRange(env, -IloInfinity, Total, 1, cName));
            GrandSubProb.add(ActiveCCSubP_UB[ActiveCCSubP_UB.getSize() - 1]);
            ExprArcsVtx.end();
            AllccVars.end();
            WhichCCs.end();
            Total.end();
        }
    }
    //Add new columns
    NewCycleTPH = IloArray<IloNumColumn>(env, PosNewCycles.size());
    for (int i = 0; i < PosNewCycles.size(); i++){
        NewCycleTPH[i] = IloNumColumn (env);
        NewCycleTPH[i] += ActiveCCSubP_UB[PosNewCycles[i]](1);
        NewCycleTPH[i] += ActiveCCSubP_LB[PosNewCycles[i]](1);
        NewCycleTPH[i] += ObjGrandSubP(GrandSubSolSet[PosNewCycles[i]].get_w());
    }
    //Generate TheOneCC
    map<int,vector<int>>::iterator it0;
    for (int i = 0; i < Pairs; i++){
        it0 = CycleNodeGSP.find(i);
        if (it0 != CycleNodeGSP.end()){
            for (int j = 0; j < CycleNodeGSP[i].size(); j++){
                int pos = FindPosVector(PosNewCycles, CycleNodeGSP[i][j]);
                if (pos != -1){
                    NewCycleTPH[pos] += TheOneCC[i](1);
                }
            }
        }
    }
    //Insert new CCs in ActiveCCSubP_LB
    map<int,int>::iterator it;
    for (int j = 0; j < GrandSubSolSet.size() - PosNewCycles.size(); j++){
        map<int,int>WhichCCs;
        for (int h = 0; h < GrandSubSolSet[j].get_cc().size(); h++){
            ccVarsOneVtx (j, GrandSubSolSet[j].get_cc()[h], WhichCCs);
        }
        for (it = WhichCCs.begin(); it != WhichCCs.end(); it++){
            for (int h = 0; h < PosNewCycles.size(); h++){
                if (it->first == PosNewCycles[h]){
                    NewCycleTPH[h] += ActiveCCSubP_LB[j](1);
                }
            }
        }
    }
    for (int i = 0; i < PosNewCycles.size(); i++){
        r.add(IloNumVar(NewCycleTPH[i]));
        SetName1Index(r[r.getSize() - 1], "r", PosNewCycles[i]);
    }
    //
}
void Problem::AddRow_MakeOneFail(){
    IloExpr ArcVerSum (env, 0);
    ArcVerSum = GenerateMakeOneFailConstraint(int(RepSolCounter));
    string name = "MakeOneFail." + to_string(RepSolCounter);
    const char* cName = name.c_str();
    ActiveGrandSubSol.add(IloRange(env, 1, ArcVerSum, IloInfinity, cName));
    GrandSubProb.add(ActiveGrandSubSol[ActiveGrandSubSol.getSize() - 1]);
}
void Problem::FillRepSol_AND_CycleNodeGSP(map<int,vector<int>>& CycleNodeGSP){
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        RepairedListCCs.push_back(Cycles(GrandSubSolSet[i].get_cc(), GrandSubSolSet[i].get_w()));
        for (int j = 0; j < RepairedListCCs.back().get_c().size(); j++){
            CycleNodeGSP[RepairedListCCs.back().get_c()[j]].push_back(i);
        }
    }
}
void Problem::ccVarsOneVtx (int cycleIndx, int vtx, map<int,int>&WhichCCs){
    IloExpr Expr (env, 0);
    map<int,int>::iterator it;
    
    for (int i = 0; i < CycleNodeGSP[vtx].size(); i++){
        if (CycleNodeGSP[vtx][i] != cycleIndx){
            it = WhichCCs.find(CycleNodeGSP[vtx][i]);
            if (it == WhichCCs.end()) WhichCCs[CycleNodeGSP[vtx][i]] = -1;
        }
    }
}
IloExpr Problem::GenerateAllCCsVars(map<int,int>&WhichCCs){
    IloExpr Expr (env,0);
    map<int,int>::iterator it;
    for (it = WhichCCs.begin(); it != WhichCCs.end(); it++){
        Expr += r[it->first];
    }
    return Expr;
}
IloExpr Problem::GenerateAllArcsVertices(int cycleIndx){
    IloExpr Expr (env, 0);

    for (int j = 0; j < GrandSubSolSet[cycleIndx].get_cc().size(); j++){
    //Find the j position in the Adjancecy matrix
        int pos = -1;
        if (j <= GrandSubSolSet[cycleIndx].get_cc().size() - 2){
            pos = FindPosArray(AdjacencyList[GrandSubSolSet[cycleIndx].get_cc()[j]], GrandSubSolSet[cycleIndx].get_cc()[j + 1] + 1);
        }
        else{
            pos = FindPosArray(AdjacencyList[GrandSubSolSet[cycleIndx].get_cc()[j]], GrandSubSolSet[cycleIndx].get_cc()[0] + 1);
        }
        Expr +=vertex[GrandSubSolSet[cycleIndx].get_cc()[j]];
        Expr += arc[GrandSubSolSet[cycleIndx].get_cc()[j]][pos];
    }
    //cout << Expr << endl;
    return Expr;
}
IloExpr Problem::GenerateMakeOneFailConstraint(int iteration){
    IloExpr ArcVerSum (env, 0);
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        //Check whether this cycle/chain was used in this iteration
        if (FindPosVector(GrandSubSolSet[i].get_itev(), iteration) != -1){
            for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
                //Find the j position in the Adjancecy matrix
                int pos = -1;
                if (j <= GrandSubSolSet[i].get_cc().size() - 2){
                    pos = FindPosArray(AdjacencyList[GrandSubSolSet[i].get_cc()[j]], GrandSubSolSet[i].get_cc()[j + 1] + 1);
                }
                else{
                    pos = FindPosArray(AdjacencyList[GrandSubSolSet[i].get_cc()[j]], GrandSubSolSet[i].get_cc()[0] + 1);
                }
                ArcVerSum +=vertex[GrandSubSolSet[i].get_cc()[j]];
                ArcVerSum += arc[GrandSubSolSet[i].get_cc()[j]][pos];
            }
            //cout << CCSum << endl;
            //name = "ActiveCycle." + to_string(i);
            //cName = name.c_str();
            //ActiveCCSubP.add(IloRange(env, 0, r[i] - 1 + CCSum, IloInfinity, cName));
        }
    }
    return ArcVerSum;
}
void Problem::CheckNewIncludedVerticesGSP(vector<IndexGrandSubSol>&Sol, vector<Cycles>&RepairedListCCs){//Modify VertexInSomeSolGSP

    //Check whether the repaired solution includes one of the previously selected (active) cycles
    for (int i = 0; i < RepairedListCCs.size(); i++){
        //check if c is in Sol already, if so add the new iteration in which is used
        int ans = checkIfCCisNew(RepairedListCCs[i].get_c(), Sol);
        if ( ans != -1){
            Sol[ans].set_ite(int(RepSolCounter));
        }
        else{//It is a new cycle/chain
            Sol.push_back(IndexGrandSubSol(RepairedListCCs[i].get_c(), RepairedListCCs[i].get_w(), int(RepSolCounter)));
            for (int j = 0; j < RepairedListCCs[i].get_c().size(); j++){
                CycleNodeGSP[RepairedListCCs[i].get_c()[j]].push_back(int(Sol.size() - 1));
            }
        }
    }
    
}
int checkIfCCisNew(vector<int>v, vector<IndexGrandSubSol>&Sol){
    
    for (int i = 0; i < Sol.size(); i++){
        int counter = 0;
        if (Sol[i].get_cc().size() == v.size() && Sol[i].get_cc()[0] <= v[0]){
            for (int j = 0; j < Sol[i].get_cc().size(); j++){
                if (Sol[i].get_cc()[j] == v[j]){
                    counter++;
                }
                else{
                    break;
                }
            }
            if (counter == v.size()) return i;
        }
    }
    return -1;//Not found
    
}
void Problem::AddNewColsConsGSP(vector<Cycles>& RepairedSol){
    //Add Bounding constraint
    IloExpr w (env,0);
    for (int i = 0; i < RepairedSol.size(); i++){
        w += RepairedSol[i].get_Many()*r[i];
    }
    string name = "BoundingC." + to_string(RepSolCounter);
    const char* cName = name.c_str();
    BoundObjGrandSubP.add(IloRange(env, 0, Beta - w, IloInfinity, cName));
    GrandSubProb.add(BoundObjGrandSubP);
    
    //Adding making-one fail
    IloExpr ArcSum (env, 0);
    IloExpr VertexSum (env, 0);
    //Add making fail one vertex/arc
    for (int i = 0; i < RepairedSol.size(); i++){
        IloExpr CCSum (env, 0);
        //IloInt nElCCSum = 0;
        for (int j = 0; j < RepairedSol[i].get_c().size(); j++){
            //Find the j position in the Adjancecy matrix
            int pos = -1;
            if (j <= RepairedSol[i].get_c().size() - 2){
                pos = FindPosArray(AdjacencyList[RepairedSol[i].get_c()[j]], RepairedSol[i].get_c()[j + 1] + 1);
            }
            else{
                pos = FindPosArray(AdjacencyList[RepairedSol[i].get_c()[j]], RepairedSol[i].get_c()[0] + 1);
            }
            // Compute for CC
            //nElCCSum+=2; //arc and vertex
            CCSum += vertex[RepairedSol[i].get_c()[j]] + arc[RepairedSol[i].get_c()[j]][pos];
            //Compute All Vertex + Arcs Sum
            VertexSum +=vertex[RepairedSol[i].get_c()[j]];
            ArcSum += arc[RepairedSol[i].get_c()[j]][pos];
        }
        //cout << CCSum << endl;
        name = "ActiveCycle." + to_string(i);
        cName = name.c_str();
        ActiveCCSubP_LB.add(IloRange(env, 0, r[i] - 1 + CCSum, IloInfinity, cName));
    }
    
    
    
}
int FindPosArray(IloNumArray array, int value){
    int pos = -1;
    for (int j = 0; j < array.getSize(); j++){
        if (array[j] == value) return j;
    }
    return pos;
}
int FindPosVector(vector<int> array, int value){
    int pos = -1;
    for (int j = 0; j < array.size(); j++){
        if (array[j] == value) return j;
    }
    return pos;
}
void Problem::SetName2Index(IloNumVar& var, const char* prefix, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
vector<Cycles> Problem::BackRecoursePolicy(IloNumArray& vertex_sol, IloNumArray2& arc_sol){
    
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        IloNumArray2 AdjaList (env, AdjacencyList.getSize());

        for (int i = 0; i < AdjacencyList.getSize(); i++){
            AdjaList[i] = IloNumArray(env);
        }
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            if (vertex_sol[GrandSubSolSet[i].get_cc()[j]] < 0.9){//If vertex did not fail
                for (int l = 0; l < AdjacencyList[GrandSubSolSet[i].get_cc()[j]].getSize(); l++){
                    int neighbour = AdjacencyList[GrandSubSolSet[i].get_cc()[j]][l] - 1;
                    int pos = FindPosVector(GrandSubSolSet[i].get_cc(), neighbour);
                    if (pos != -1) {
                        if (arc_sol[GrandSubSolSet[i].get_cc()[j]][l] < 0.9){//If arc did not fail
                            AdjaList[GrandSubSolSet[i].get_cc()[j]].add(neighbour + 1);
                        }
                    }
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
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            int origin = GrandSubSolSet[i].get_cc()[j];
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
    return CFThirdPhase(NewListCCs, CycleNodeTPH);
}
vector<Cycles> Problem::AmongPolicy(IloNumArray& vertex_sol, IloNumArray2& arc_sol){
    
    IloNumArray2 AdjaList (env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env);
    }
    
    vector<int>LookFor;
    
    for (int i = 0; i < GrandSubSolSet.size(); i++){
        for (int j = 0; j < GrandSubSolSet[i].get_cc().size(); j++){
            if (vertex_sol[GrandSubSolSet[i].get_cc()[j]] < 0.9){//If vertex did not fail
                LookFor.push_back(GrandSubSolSet[i].get_cc()[j]);
            }
        }
    }
    
    //Fill in Adjacency List
    for (int i = 0; i < LookFor.size(); i++){
        //cout << endl << LookFor[i] << "\t" << ":";
        for (int j = 0; j < AdjacencyList[LookFor[i]].getSize(); j++){
            int pos = FindPosVector(LookFor, AdjacencyList[LookFor[i]][j] - 1);
            if (pos != -1){
                if (arc_sol[LookFor[i]][j] < 0.9){//If arc did not fail
                    AdjaList[LookFor[i]].add(AdjacencyList[LookFor[i]][j]);
                    //cout << AdjacencyList[LookFor[i]][j] - 1 << "\t";
                }
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
            for (int l = 0; l < NewList[k].get_c().size(); l++){
                CycleNodeTPH[NewList[k].get_c()[l]].push_back(int(NewListCCs.size() - 1));//CycleNodeMap for the third phase
            }
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    
    return CFThirdPhase(NewListCCs, CycleNodeTPH);
}
vector<Cycles> Problem::AllPolicy(IloNumArray& vertex_sol, IloNumArray2& arc_sol){
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
    return CFThirdPhase(NewListCCs, CycleNodeTPH);
}
vector<Cycles> Problem::CFThirdPhase(vector<Cycles>Input, map<int,vector<int>>&CycleNodeMap){
    vector<Cycles>SolThirdPH;

    IloEnv env2;
    ThirdPH = IloModel(env2);
    cplexThirdPH = IloCplex(ThirdPH);
    cplexThirdPH.setParam(IloCplex::Param::TimeLimit, 1800);
    cplexThirdPH.setParam(IloCplex::Param::Threads, 1);
    
    //Create variables
    IloNumVarArray y(env2, Input.size(), 0, 1, ILOINT);
    for (int i = 0; i < Input.size(); i++){
       SetName1Index(y[i], "y", i);
       //cout << y[i].getName() << endl;
    }

    //Add constraints
    for (int i = 0; i < Pairs; i++){
       IloExpr cycleTPH (env2,0);
       auto pos = CycleNodeMap.find(i);
       if (pos != CycleNodeMap.end()){
           for (int j = 0; j < CycleNodeMap[i].size(); j++){
               cycleTPH+= y[CycleNodeMap[i][j]];
           }
           //Add constraint: node only in one cycle
           ThirdPH.add(cycleTPH <= 1);
       }
       cycleTPH.end();
    }
    //Objective value
    IloExpr ObjTPH(env2, 0);
    for (int j = 0; j < Input.size(); j++){
       ObjTPH += Input[j].get_w()*y[j];
    }
    ThirdPH.add(IloMaximize(env2, ObjTPH));
    //cplexThirdPH.exportModel("ThirdPhase.lp");
    
    cplexThirdPH.solve();
    IloNumArray y_sol(env2, Input.size());
    cplexThirdPH.getValues(y_sol,y);
    
    env2.out() << "TPH Objective: " << cplexThirdPH.getObjValue() << endl;

    for (int j = 0; j < y_sol.getSize(); j++){
        if (y_sol[j] > 0.9){
            SolThirdPH.push_back(Input[j]);
            SolThirdPH.back().set_Many(int(Input[j].get_c().size()));
        }
    }
    ThirdPH.end();
    cplexThirdPH.end();
    env2.end();
    
    return SolThirdPH;
}
