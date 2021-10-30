//
//  KEP_Deterministic.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-16.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandProblem.hpp"
void Problem::KEP_CycleFormulation(){
    //Create model
        GrandProb = IloModel (env);
        cplexGrandP = IloCplex(GrandProb);
        cplexGrandP.setParam(IloCplex::Param::TimeLimit, 1800);
        cplexGrandP.setParam(IloCplex::Param::Threads, 1);
        IloNum time1 = cplexGrandP.getCplexTime();

        //Find Cylces
        MainCycleFinder();

        //Create variables
        z = IloNumVarArray(env, ListCycles.size(), 0, 1, ILOINT);
        for (int i = 0; i < ListCycles.size(); i++){
           SetName1Index(z[i], "z", i);
           //cout << z[i].getName() << endl;
        }
        eta = IloNumVar(env);
        SetName1Index(eta, "eta", 0);
    
//        for (int i = 0; i < ListCycles.size(); i++){
//           SetName1Index(z[i], "z", i);
//           //cout << z[i].getName() << endl;
//        }

        //Add constraints
        for (int i = 0; i < Pairs; i++){
           IloExpr cycle (env,0);
           auto pos = CycleNode.find(i);
           if (pos != CycleNode.end()){
               for (int j = 0; j < CycleNode[i].size(); j++){
                   cycle+= z[CycleNode[i][j]];
               }
           }
           //Add constraint: node only in one cycle
           GrandProb.add(cycle <= 1);
           cycle.end();
        }

        //Objective value
        IloExpr Obj(env, 0);
        for (int j = 0; j < ListCycles.size(); j++){
           Obj += ListCycles[j].get_w()*z[j];
        }
        GrandProb.add(eta <= Obj);
        GrandProb.add(IloMaximize(env, eta));
    
        //Cutting plane
        ActiveGrandSubSol = IloRangeArray(env);
        BoundObjective = IloRangeArray(env) ;
        IloBool NewScenario = true;
    
        while(NewScenario == true){
            Iteration++;
            cplexGrandP.solve();
            IloNum time2 = cplexGrandP.getCplexTime();

           if (cplexGrandP.getStatus() == IloAlgorithm::Infeasible) {
               env.out() << "No solution" << endl;
           }
           else {
               CFObj = cplexGrandP.getObjValue();
               BestObj = cplexGrandP.getBestObjValue();
               Gap = cplexGrandP.getMIPRelativeGap();
               NCycles = ListCycles.size();
               SolTime = time2 - time1;
               
               env.out() << "Objective: " << CFObj << endl;
               
               IloNumArray z_sol(env, ListCycles.size());
               cplexGrandP.getValues(z_sol,z);
               
               GrandProbSol.clear();
               for (int f = 0; f < z_sol.getSize(); f++){
                   if (z_sol[f] > 0.9){
                       GrandProbSol.push_back(IndexGrandSubSol(ListCycles[f].get_c(), ListCycles[f].get_w()));
                       GrandProbSol.back().set_ite(1);
                   }
               }
               //GrandSubProbMaster();
               cout << "hi";
           }
        }
    
}
void Problem::AddNewColsConsGP(){
    //Find cycle index
    vector<int> ref;
    ref.push_back(1);
    ref.push_back(37);
    ref.push_back(35);
    
    auto pos = CycleNode.find(1);
    int index = -1;
    if (pos != CycleNode.end()){
        for (int i = 0; i < CycleNode[1].size(); i++){
            int count = 0;
            for (int j = 0; j < ListCycles[CycleNode[1][i]].get_c().size(); j++){
                for (int r = 0; r < ref.size(); r++){
                    //cout << ListCycles[CycleNode[1][i]].get_c()[j] << "\t" << CycleNode[1][i] << endl;
                    if (ListCycles[CycleNode[1][i]].get_c()[j] == ref[r]){
                        count++;
                        break;
                    }
                }
                if (count == ref.size()){
                    index = CycleNode[1][i];
                    break;
                }
            }
            if (count == ref.size()) break;
        }
    }

    //Add new constraints
    GrandSubOptSol.push_back(7);
    if (GrandSubOptSol.size() > 0){
        for (int i = 0; i < GrandSubOptSol.size(); i++){
            string name = "RecourseSol_PlannedSol." + to_string(i) + "." + to_string(Iteration);
            const char* cName = name.c_str();
            ActiveGrandSubSol.add(IloRange(env, -IloInfinity, -z[index], 0, cName));
        }
        GrandProb.add(ActiveGrandSubSol);
        //Add new columns
        IloExpr newTotal (env,0);
        for (int i = 0; i < GrandSubOptSol.size(); i++){
            cgrandpsol = IloNumColumn (env);
            cgrandpsol += ActiveGrandSubSol[i](1);
            IloNumVar newcc(cgrandpsol);
            SetName1Index(newcc, "x", i);
            newTotal-= ListCycles[GrandSubOptSol[i]].get_w()*newcc;
        }
        newTotal+= eta;
        string name = "BoundGrandP." + to_string(Iteration);
        const char* cName = name.c_str();
        BoundObjective.add(IloRange(env, newTotal, 0, cName));
        cout << BoundObjective[0];
        GrandProb.add(BoundObjective);
    }
    cplexGrandP.exportModel("GrandProb.lp");
    cplexGrandP.solve();
    CFObj = cplexGrandP.getObjValue();
    env.out() << "Objective: " << CFObj << endl;
    
}
void Problem::SetName1Index(IloNumVar& var, const char* prefix, IloInt i){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i);
    const char* varName = name.c_str();
    var.setName(varName);
}
void Problem::MainCycleFinder(){
    IloModel m(env);
    IloCplex cp(m);
    
    //Pick a source node
    int origin;
    
    //for (int u = 0; u < comp.size(); u++){
        //Copy Adjacency List
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    for (int l = 0; l < AdjacencyList.getSize(); l++){
        AdjaList[l] = IloNumArray(env);
        for (int i = 0; i < AdjacencyList[l].getSize(); i++){
            AdjaList[l].add(AdjacencyList[l][i]);
        }
        //cout << AdjaList[l].getSize() << " " << AdjacencyList[l].getSize() << endl;
    }
    bool terminate = false;
    origin = -1;
    while (terminate == false){
        //Set origin
        int u = origin + 1;
        int prev_origin = origin;
        for (; u < Pairs; u++){
            if (AdjaList[u].getSize() > 0){
                origin = u;//subtracted 1, so 100 is 99
                break;
            }
        }
        if (origin == prev_origin){
            break;
        }
        
        //Find Cycles
        int ThisMany = int(ListCycles.size());
        vector<Cycles> ListC;
        ListC = SubCycleFinder(env, AdjaList, origin);
        for (int i = 0; i < ListC.size(); i++){
            ListCycles.push_back(ListC[i]);
            for (int y = 0; y < ListC[i].get_c().size(); y++){
                CycleNode[ListC[i].get_c()[y]].push_back(ListCycles.size() - 1);
            }
        }
        if (ListCycles.size() == ThisMany) terminate = true;
        
        //Remove node
        for (int l = 0; l < AdjaList.getSize(); l++){
            if (l == origin){
                //cout << AdjaList[l].getSize() << endl;
                AdjaList[l] = IloNumArray(env,0);
                //cout << AdjaList[l].getSize() << endl;
            }
            else{
                for (int i = 0; i < AdjaList[l].getSize(); i++){
                    if (AdjaList[l][i] == origin + 1){
                        AdjaList[l].remove(i);
                    }
                }
            }
        }
        
    }
    
    m.end();
    cp.end();
    
}
vector<Cycles> Problem::SubCycleFinder (IloEnv env, IloNumArray2 AdjaList, IloInt origin){
    vector<Cycles> ListC;
    vector<int> nodesBestCycle;
    IloNumArray stepinSU (env, AdjaList.getSize());
    for (int c = 0; c < stepinSU.getSize(); c++) stepinSU[c] = 0;
    vector<int> Stack;
    Stack.push_back(origin);
    IloInt whoi = origin;
    stepinSU[origin] = 0;
    IloBool VLpresent = true; //because we start with the origin
    IloBool cont = true;
    IloInt prevWhoi = 1000000;
    IloNum BigW = 0;
    IloNum weight = 0;
    IloBool NewCycle = false;
    if (AdjaList[whoi].getSize() == 0){
        cont = false;
    }
    //cout << endl << dicX[sp][this->vL[sp] - 1] << " " << this->vL[sp];
    while (cont == true){
        //cout << AdjaList[whoi][stepinSU[whoi]] << endl;
        if (IsxinStack(AdjaList[whoi][stepinSU[whoi]] - 1, Stack) == true || Stack.size() >= CycleLength){
            if (AdjaList[whoi][stepinSU[whoi]] == (origin + 1) && VLpresent == true && NewCycle == false){
                // they form a cycle together!
                weight = PathWeight(Stack);
                ListC.push_back(Cycles(Stack, weight));
//                for (int y = 0; y < Stack.size(); y++){
//                    CycleNode[Stack[y]].push_back(ListCycles.size() - 1);
//                }
                if (weight > BigW){
                    BigW = weight;
                    nodesBestCycle = Stack;
                }
                //Backtrack
                NewCycle = true;
            }
            else{
                NewCycle = false;
                while(true){
                    //Backtrack
                    stepinSU[whoi]++;
                    if (stepinSU[whoi] >= AdjaList[whoi].getSize()){
                        if (Stack.back() == origin){ VLpresent = false;}
                        stepinSU[whoi] = 0;
                        Stack.pop_back();
                        if (Stack.size() == 0){
                            cont  = false;
                            break;
                        }else{
                            whoi = Stack.back();
                        }
                    }
                    else{
                        break;
                    }
                }
            }
        }
        else{
            if (AdjaList[whoi][stepinSU[whoi]] == origin + 1) VLpresent = true;
            Stack.push_back(AdjaList[whoi][stepinSU[whoi]] - 1);
            prevWhoi = whoi;
            whoi = AdjaList[whoi][stepinSU[whoi]] - 1;
            if (AdjaList[whoi].getSize() == 0){
                //bactrack
                if (Stack.back() == origin){ VLpresent = false;}
                stepinSU[whoi] = 0;
                Stack.pop_back();
                whoi = Stack.back();
                stepinSU[whoi]++;
                while(true){
                    if (stepinSU[whoi] >= AdjaList[whoi].getSize()){
                        if (Stack.back() == origin){VLpresent = false;}
                        stepinSU[whoi] = 0;
                        Stack.pop_back();
                        if (Stack.size() == 0){
                            cont  = false;
                            //return BigW;
                        }else{
                            whoi = Stack.back();
                            stepinSU[whoi]++;
                        }
                    }
                    else{
                        break;
                    }
                }
            }
        }
    }
    return ListC;
}
IloNum Problem::PathWeight (vector<int>& Stack){
    IloNum weight = 0;
    //arc weights
    for (int l = 0; l < Stack.size(); l++){
        if (l >= 1){
            for (int i = 0; i < PredList[Stack[l]].size(); i++){
                if (PredList[Stack[l]][i] == Stack[l - 1]){
                    weight += Weights[make_pair(Stack[l - 1], Stack[l])];
                    break;
                }
            }
        }
        else{
            for (int i = 0; i < PredList[Stack[l]].size(); i++){
                if (PredList[Stack[l]][i] == Stack[Stack.size() - 1]){
                    weight += Weights[make_pair(Stack[Stack.size() - 1], Stack[l])];
                    //weight += PredList[Stack[l]][i].second;
                    break;
                }
            }
        }
    }
    return weight;
}
IloBool Problem::IsxinStack (IloInt test, vector<int>& xinTrial){
    for (int i = 0; i < xinTrial.size(); i++){
        if (test == xinTrial[i])    return true;
    }
    return false;
}
