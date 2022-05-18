//
//  KEP_Deterministic.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-16.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandProblem.hpp"
#include <random>
random_device rd;
mt19937 gen(rd());


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
                CycleNode[ListC[i].get_c()[y]].push_back(int(ListCycles.size()) - 1);
                if (y < ListC[i].get_c().size() - 1){
                    CycleArcs[make_pair(ListC[i].get_c()[y], ListC[i].get_c()[y + 1])].push_back(int(ListCycles.size()) - 1);
                }
                else{
                    CycleArcs[make_pair(ListC[i].get_c()[y], ListC[i].get_c()[0])].push_back(int(ListCycles.size()) - 1);
                }
            }
        }
        //if (ListCycles.size() == ThisMany) terminate = true;
        
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



//Creación de variables
NumVar4D Create4DBin (IloEnv& env, IloNumArray2 Adja, int L, int U, string name){
    NumVar4D Aux(env, Adja.getSize());
    for (int i = 0; i < Aux.getSize(); i++) {
        Aux[i] = NumVar3D(env, Adja[i].getSize());
        for (int j = 0; j < Aux[i].getSize(); j++) {
            Aux[i][j] = NumVar2D(env, L);
            for (int l = 0; l < L; l++) {
                Aux[i][j][l] = IloNumVarArray(env, U, 0, 1, ILOINT);
                for (int u = 0 ; u < U; u++) {
                    string auxName = name + "_" + to_string(i) + "_" + to_string(int(Adja[i][j] - 1)) + "_" + to_string(l) + "_" + to_string(u);
                    Aux[i][j][l][u].setName(auxName.c_str());
                }
            }
        }
    }
    return Aux;
}
void SetUB4DBin (IloEnv& env, NumVar4D& vars, vector<int>& distNDD, int P, int L, int u){
    for (int i = 0; i < vars.getSize(); i++) {
        if (i < P){//Only patient-donor pairs
            for (int j = 0; j < vars[i].getSize(); j++) {
                for (int l = 0; l < L; l++) {
                    if (distNDD[i] + 1 > L || l == 0){//
                        vars[i][j][l][u].setUB(0);
                    }
                }
            }
        } else {
            for (int j = 0; j < vars[i].getSize(); j++) {
                for (int l = 1; l < L; l++) {//Arc (i,j) where i is an altruist, can only be at position 0 of a chain
                    vars[i][j][l][u].setUB(0);
                }
            }
        }
    }
}
NumVar3D Create3DBin (IloEnv& env, IloNumArray2 Adja, int L, string name){
    NumVar3D Aux(env, Adja.getSize());
    for (int i = 0; i < Adja.getSize(); i++) {
        Aux[i] = NumVar2D(env, Adja[i].getSize());
        for (int j = 0; j < Adja[i].getSize(); j++) {
            Aux[i][j] = IloNumVarArray(env, L, 0, 1, ILOINT);
            for (int l = 0; l < L; l++) {
                string auxName = name + "_" + to_string(i + 1) + "_" + to_string(int(Adja[i][j])) + "_" + to_string(l);
                Aux[i][j][l].setName(auxName.c_str());
            }
        }
    }
    
    return Aux;
}
void SetUB3DBin (IloEnv& env, NumVar3D& vars, vector<int>& distNDD, int P, int L){
    for (int i = 0; i < vars.getSize(); i++) {
        if (i < P){//Only patient-donor pairs
            for (int j = 0; j < vars[i].getSize(); j++) {
                for (int l = 0; l < L; l++) {
                    if (distNDD[i] + 1 > L || l == 0){
                        vars[i][j][l].setUB(0);
                    }
                }
            }
        } else {
            for (int j = 0; j < vars[i].getSize(); j++) {
                for (int l = 1; l < L; l++) {//Arc (i,j) where i is an altruist, can only be at position 0 of a chain
                    vars[i][j][l].setUB(0);
                }
            }
        }
    }
}
NumVar2D Create2DBin (IloEnv& env, int I, int J, string name){
    NumVar2D Aux(env, I);
    for (int i = 0; i < I; i++) {
        Aux[i] = IloNumVarArray(env, J, 0, 1, ILOINT);
        for (int j = 0; j < J; j++) {
            string auxName = name + "_" + to_string(i) + "_" + to_string(j);
            Aux[i][j].setName(auxName.c_str());
        }
    }
    return Aux;
}
IloNumVarArray Create1DBin (IloEnv& env, int I, string name){
    IloNumVarArray Aux = IloNumVarArray(env, I, 0, 1, ILOINT);
    for (int i = 0; i < I; i++) {
        string auxName = name + "[" + to_string(i) + ']';
        Aux[i].setName(auxName.c_str());
    }
    return Aux;
}

//Parameters creation
void CreateU_uj(int maxFailureV, int Nodes, int maxFailureA, vector<map<pair<int,int>, bool>>&scenarios, IloNumArray2 AdjaList){
    scenarios.push_back(map<pair<int,int>, bool>());
    uniform_int_distribution<int> distri(0, Nodes - 1);
    vector<int>sceV;
    int cuantos = 0;
    if (maxFailureV > 0){
        while (cuantos < maxFailureV){
            int cual = Nodes + rand()%Nodes; // Artificial scenario
            //int cual = distri(gen);
            bool isin = false;
            for (int i = 0; i < sceV.size(); i++){
                if (sceV[i] == cual) isin = true;
            }
            if (isin == false){
                scenarios.back()[make_pair(-1,cual)] = true;
                cuantos++;
            }
        }
    }

    cuantos = 0;
    if (maxFailureA > 0){
        while (cuantos < maxFailureA){
            int cual = Nodes + rand()%Nodes; // Artificial scenario
            bool isin = false;
            for (int i = 0; i < sceV.size(); i++){
                if (sceV[i] == cual) isin = true;
            }
            if (isin == false){
                //uniform_int_distribution<int> distri2(0, int(AdjaList[cual].getSize() -1));
                int cual2 = Nodes + rand()%Nodes; // Artificial scenario
                //int cual2 = distri2(gen);
                scenarios.back()[make_pair(cual, cual2)] = true;
                cuantos++;
            }
        }
    }
}

//Creación de restricciones
IloRangeArray CreateCons7b(IloEnv& env, IloNumVar& Z, NumVar2D& Y_ju, int Ulast, int Pairs, string name){
    IloRangeArray Cons(env);
    IloExpr expr (env, 0);
    expr += Z;
    for (int j = 0; j < Pairs; j++){
        expr -= Y_ju[j][Ulast];
    }
    string name2 = name + '[' + to_string(Ulast) + ']';
    Cons.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
    return Cons;
}
IloRangeArray CreateCons7c(IloEnv& env, NumVar2D& Y_ju, IloNumVarArray& X_c, NumVar3D& E_ijl, int Ulast, int Pairs, map<int,vector<int>> CycleNode, map<int,vector<pair<int,int>>> PredMap, string name){
    IloRangeArray Cons(env);
    for (int j = 0; j < Pairs; j++) {   // para tood j
        IloExpr expr (env, 0);
        expr+= Y_ju[j][Ulast];
        for (int c = 0; c < CycleNode[j].size(); c++) { // - suma C_kj
            expr-= X_c[CycleNode[j][c]];  // pongo el c correspondiente de CycleNode de un j
        }
        for (int i = 0; i < PredMap[j].size(); i++) {
            int i2 = PredMap[j][i].first;
            int j2 = PredMap[j][i].second;
            expr-= IloSum(E_ijl[i2][j2]); // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
        }
        string name2 = name + '_' + to_string(Ulast) + '_' + to_string(j);
        Cons.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
    }
    return Cons;
}
IloRangeArray CreateCon7d(IloEnv& env, NumVar2D& Y_ju, NumVar2D& X_cu, NumVar4D& E_ijlu, int Ulast, int Pairs, map<int,vector<int>> CycleNode, int Lmax, map<int,vector<pair<int,int>>> PredMap, string name){
    IloRangeArray Cons(env);
    for (int j = 0; j < Pairs; j++) { // para todo j
        IloExpr expr (env, 0);
        expr+= Y_ju[j][Ulast];
        for (int c = 0; c < CycleNode[j].size(); c++) { // - suma C_kj
            expr-= X_cu[CycleNode[j][c]][Ulast];  // pongo el c correspondiente de CycleNode de un j
        }
        for (int i = 0; i < PredMap[j].size(); i++) {
            for (int l = 0; l < Lmax; l++){
                int i2 = PredMap[j][i].first;
                int j2 = PredMap[j][i].second;
                expr-= E_ijlu[i2][j2][l][Ulast]; // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
            }
        }
        string name2 = name + '[' + to_string(Ulast) + ']' + '[' + to_string(j) + ']';
        Cons.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
        expr.end();
    }
    return Cons;
}
IloRangeArray CreateCon7e(IloEnv& env, IloNumVarArray& var2D, NumVar3D& var3D, int Pairs, map<int,vector<int>> CycleN, map<int,vector<pair<int,int>>> PredMap, string name){
    IloRangeArray Cons(env);
    for (int j = 0; j < Pairs; j++) { // para todo j
        IloExpr expr (env, 0);
        for (int c = 0; c < CycleN[j].size(); c++) {
            expr+= var2D[CycleN[j][c]];
        }
        for (int i = 0; i < PredMap[j].size(); i++) {
            int i2 = PredMap[j][i].first;
            int j2 = PredMap[j][i].second;
            expr+= IloSum(var3D[i2][j2]);
        }
        string name2 = name + '[' + to_string(j) + ']';
        Cons.add(IloRange(env, -IloInfinity, expr, 1, name2.c_str()));
        expr.end();
    }
    return Cons;
}
IloRangeArray CreateCon7f(IloEnv& env, NumVar3D& var3D, int P, IloNumArray2 Adja, string name){
    IloRangeArray Cons(env);
    for (int j = P; j < Adja.getSize(); j++) { // para todo j
        IloExpr expr (env, 0);
        for (int i = 0; i < Adja[j].getSize(); i++) {
            expr+= var3D[j][i][0];
        }
        string name2 = name + '[' + to_string(j) + ']';
        Cons.add(IloRange(env, -IloInfinity, expr, 1, name2.c_str()));
        expr.end();
    }
    return Cons;
}
IloRangeArray CreateCon7g(IloEnv& env, NumVar2D& X_cu, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, int Pairs, map<int,vector<int>> CycleNode, map<int, vector<pair<int,int>>> PredMap, int Lmax, string name){
    IloRangeArray Cons(env);
    for (int j = 0; j < Pairs; j++) { // para todo j
        IloExpr expr (env, 0);
        auto it = scenarios[Ulast].find(make_pair(-1, j));
        if (it != scenarios[Ulast].end()){
            for (int c = 0; c < CycleNode[j].size(); c++) {
                X_cu[CycleNode[j][c]][Ulast].setUB(0);
            }
            for (int l = 0; l < Lmax; l++) {
                for (int i = 0; i < PredMap[j].size(); i++) {
                    int i2 = PredMap[j][i].first;
                    int j2 = PredMap[j][i].second;
                    //cout << E_ijlu[i2][j2][l][Ulast].getName();
                    E_ijlu[i2][j2][l][Ulast].setUB(0); // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
                }
            }
        }
        else{
            for (int c = 0; c < CycleNode[j].size(); c++) {
                expr+= X_cu[CycleNode[j][c]][Ulast];
            }
            for (int l = 0; l < Lmax; l++) {
                for (int i = 0; i < PredMap[j].size(); i++) {
                    int i2 = PredMap[j][i].first;
                    int j2 = PredMap[j][i].second;
                    expr+= E_ijlu[i2][j2][l][Ulast]; // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
                }
            }
            string name2 = name + '[' + to_string(Ulast) + ']' + '[' + to_string(j) + ']';
            Cons.add(IloRange(env, -IloInfinity, expr, 1, name2.c_str()));
        }
        
    }
    return Cons;
}
IloRangeArray CreateCon7h(IloEnv& env, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, int P, IloNumArray2 Adja, string name){
    IloRangeArray Cons(env);

    for (int j = P; j < Adja.getSize(); j++) { // para todo j
        IloExpr expr (env, 0);
        auto it = scenarios[Ulast].find(make_pair(-1, j));
        if (it != scenarios[Ulast].end()){
            for (int i = 0; i < Adja[j].getSize(); i++) {
                E_ijlu[j][i][0][Ulast].setUB(0);
            }
        }
        else{
            for (int i = 0; i < Adja[j].getSize(); i++) {
                expr+= E_ijlu[j][i][0][Ulast];
            }
            string name2 = name + '[' + to_string(Ulast) + ']' + '[' + to_string(j) + ']';
            Cons.add(IloRange(env, -IloInfinity, expr, 1, name2.c_str()));
        }
    }
    return Cons;
}
IloRangeArray CreateCon7i(IloEnv& env, NumVar3D& var3D, int Pairs, int Lmax, map<int, vector<pair<int,int>>> PredMap, IloNumArray2 Adja, string name){
    IloRangeArray Cons(env);
    for (int j = 0; j < Pairs; j++) { // para todo j
        for (int l = 1; l < Lmax; l++) { // para todo l, pero L \ {1}
            IloExpr expr (env, 0);
            for (int i = 0; i < PredMap[j].size(); i++) {
                int i2 = PredMap[j][i].first;
                int j2 = PredMap[j][i].second;
                expr-= var3D[i2][j2][l - 1]; // (i,pos en la que está 8) {(4,8)  (7,8)  (10,8)}  size 3    se pone pos 0,5,10
            }
            for (int i = 0; i < Adja[j].getSize(); i++) {
                expr+= var3D[j][i][l];
            }
            //cout << expr << endl;
            string name2 = name + '[' + to_string(j) + ']' + '[' + to_string(l) + ']';
            Cons.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
        }
    }
    return Cons;
}
IloRangeArray CreateCon7j(IloEnv& env, NumVar4D& E_ijlu, int Ulast, int Pairs, int Lmax, map<int, vector<pair<int,int>>> PredMap, IloNumArray2 Adja, string name){
    IloRangeArray Cons(env);

    for (int j = 0; j < Pairs; j++) { // para todo j
        for (int l = 1; l < Lmax; l++) { // para todo l
            IloExpr expr (env, 0);
            for (int i = 0; i < PredMap[j].size(); i++) {
                int i2 = PredMap[j][i].first;
                int j2 = PredMap[j][i].second;
                expr-= E_ijlu[i2][j2][l - 1][Ulast];
            }
            for (int i = 0; i < Adja[j].getSize(); i++) {
                expr+= E_ijlu[j][i][l][Ulast];
            }
            //cout << expr << endl;
            string name2 = name + '_' + to_string(Ulast) + '_' + to_string(j) + '_' + to_string(l);
            Cons.add(IloRange(env, -IloInfinity, expr, 0, name2.c_str()));
        }
    }
    return Cons;
}
void CreateCon7k(IloEnv& env, NumVar2D& X_cu, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, IloNumArray2 AdjaList, map<pair<int,int>,vector<int>> CycleArcs, map<int, vector<pair<int,int>>> PredMap, int Lmax, string name){
    
    
    for (auto it = scenarios[Ulast].begin(); it != scenarios[Ulast].end(); it++) { // para todo j
        if (it->first.first != -1){
            for (int c = 0; c < CycleArcs[it->first].size(); c++) {
                X_cu[CycleArcs[it->first][c]][Ulast].setUB(0);
            }
            // para todo j
            for (int j = 0; j < AdjaList[it->first.first].getSize(); j++) {
                if (AdjaList[it->first.first][j] == it->first.second + 1){
                    for (int l = 0; l < Lmax; l++) {
                        E_ijlu[it->first.first][j][l][Ulast].setUB(0);
                    }
                    break;
                }
            }
        }
    }
    
}
//Retrieve cycles and chains
vector<vector<int>> GetChainsFrom1stStageSol(IloNumArray2 AdjacencyList,IloNumArray3 ysol, int Pairs, int ChainLength){
    vector<vector<int>>vChains;
    int l;
    
    map<int,bool> vertexList;
    int itgoes = Pairs;
    for (int i = itgoes; i < ysol.getSize(); i++){
        l = 0;
        vChains.push_back(vector<int>());
        vChains.back().push_back(i);
        for (int j = 0; j < ysol[i].getSize(); j++){
            if (ysol[i][j][l] > 0.9){
                vChains.back().push_back(AdjacencyList[i][j] - 1);
                //Complete chain
                i = AdjacencyList[i][j] - 1;
                j = -1;
                l++;
                if (l == ChainLength){
                    break;
                }
            }
        }
        if (vChains.back().size() == 1) vChains.erase(vChains.end() - 1);
        itgoes++;
        i = itgoes - 1;
    }
    
    
    return vChains;
}
///////////////////////////////////////////////////
//////////////////// Robust_KEP/////////////
//////////////////////////////////////////////////

void Problem::ROBUST_KEP(){
    tStart1stS = clock();
    //Dijkstra chains;
    vector<int> dist;
    distNDD = vector<int>(AdjacencyList.getSize(), INT_MAX);
    for (int i = int(Pairs); i < AdjacencyList.getSize(); i++){//Distance from altruists
        dist = dijkstra(AdjacencyList, i);
        for(int j = 0; j < dist.size(); j++){
            if (dist[j] < distNDD[j]) distNDD[j] = dist[j];//Find smallest distance from any altruist to a vertex
        }
    }
    
    //Create model
    RobustMod = IloModel (env);
    
    //Solver
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    cplexRobust = IloCplex(RobustMod);
    cplexRobust.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexRobust.setParam(IloCplex::Param::Threads, 1);
    cplexRobust.setOut(env.getNullStream());
    
    // Find cycles
    MainCycleFinder();
    
    //Create variables
    
    // Z
    string auxName = "Z";
    Z = IloNumVar(env, 0, IloInfinity, ILOFLOAT, auxName.c_str()); // pero no hay dominio en pdf
    
    // Num escenarios y U_uj
    int Umax = 1;
    CreateU_uj(int(MaxVertexFailures), int(Nodes), int(MaxArcFailures), scenarios, AdjacencyList);
    
    ////////////Variables///////////
    // Y_ju
    Y_ju = Create2DBin(env, int(Pairs), Umax, "Y_ju");
    
    // X_c
    X_c = Create1DBin (env, int(ListCycles.size()), "X_c");
    
    // X_cu
    X_cu = Create2DBin(env, int(ListCycles.size()), Umax, "X_cu");
    
    // E_ijl
    E_ijl = Create3DBin(env, AdjacencyList, int(ChainLength), "E_ijl");
    
    //Set UB on "E_ijl"
    SetUB3DBin (env, E_ijl, distNDD, int(Pairs), int(ChainLength));
    
    // E_ijlu
    E_ijlu = Create4DBin(env, AdjacencyList, int(ChainLength), Umax, "E_ijlu");
    
    SetUB4DBin(env, E_ijlu, distNDD, int(Pairs), int(ChainLength), 0);
    
    ////////////Restrictions///////////
    //Constraint (7b)
    IloRangeArray cons7b(env);
    cons7b = CreateCons7b(env, Z, Y_ju, 0, int(Pairs), "_7b");
    RobustMod.add(cons7b);
    
    //Constraint (7c)
    IloRangeArray cons7c (env);
    cons7c = CreateCons7c(env, Y_ju, X_c, E_ijl, 0, int(Pairs), CycleNode, PredMap, "7c");
    RobustMod.add(cons7c);
    
    //Constraint (7d)
    IloRangeArray cons7d (env);
    cons7d = CreateCon7d(env, Y_ju, X_cu, E_ijlu, 0, int(Pairs), CycleNode, int(ChainLength), PredMap, "7d");
    RobustMod.add(cons7d);
    
    //Constraint (7e)
    IloRangeArray cons7e (env);
    cons7e = CreateCon7e(env, X_c, E_ijl, int(Pairs), CycleNode, PredMap, "7e");
    RobustMod.add(cons7e);
    
    //Constraint (7f)
    IloRangeArray cons7f (env);
    cons7f = CreateCon7f(env, E_ijl, int(Pairs), AdjacencyList, "7f");
    RobustMod.add(cons7f);
    
    //Constraint (7g)
    IloRangeArray cons7g (env);
    cons7g = CreateCon7g(env, X_cu, E_ijlu, scenarios, 0, int(Pairs), CycleNode, PredMap, int(ChainLength), "7g");
    RobustMod.add(cons7g);
    
    //Constraint (7h)
    IloRangeArray cons7h (env);
    cons7h = CreateCon7h(env, E_ijlu, scenarios, 0, int(Pairs), AdjacencyList, "7h");//Set UB to 0 when U_uj = 1
    RobustMod.add(cons7h);
    
    //Constraint (7i)
    IloRangeArray cons7i (env);
    cons7i = CreateCon7i(env, E_ijl, int(Pairs), int(ChainLength), PredMap, AdjacencyList, "7i");
    RobustMod.add(cons7i);
    
    //Constraint (7j)
    IloRangeArray cons7j (env);
    cons7j = CreateCon7j(env, E_ijlu, 0, int(Pairs), int(ChainLength), PredMap, AdjacencyList, "7j");
    RobustMod.add(cons7j);
    
//    //Constraint (7k)
//    IloRangeArray cons7k (env);
//    CreateCon7k(env, X_cu, E_ijlu, scenarios, 0, AdjacencyList, CycleArcs, PredMap, int(ChainLength), "7k");
//    RobustMod.add(cons7k);

    //Objective function
    auxName = "ZFRuMax";
    IloObjective exprObj = IloObjective(env, Z, IloObjective::Maximize, auxName.c_str());
    RobustMod.add(exprObj);
    
    //cplexRobust.exportModel("RO_Model.lp");
    cplexRobust.solve();
    Ite1stStage++;
    //cout << "Status " << cplexRobust.getStatus() << endl;
    //cout << "Objetive: " << cplexRobust.getObjValue();
    
    
    //Retrieve solution
    FPMIP_Obj = cplexRobust.getObjValue();
    vector<IndexGrandSubSol>SolFirstStage;

//    cout << endl << "Cycles: " << endl;
    IloNumArray xsol(env, ListCycles.size());
    cplexRobust.getValues(xsol,X_c);
    for (int i = 0; i < xsol.getSize(); i++){
        if (xsol[i] > 0.9){
            SolFirstStage.push_back(IndexGrandSubSol(ListCycles[i].get_c(), ListCycles[i].get_c().size()));
        }
    }
    
//    cout << endl << "Chains: " << endl;
    IloNumArray3 esol(env, AdjacencyList.getSize());
    for (int i = 0; i < esol.getSize(); i++){
        esol[i] = IloNumArray2 (env, AdjacencyList[i].getSize());
        for (int j = 0; j < esol[i].getSize(); j++){
            esol[i][j] = IloNumArray(env, ChainLength);
            cplexRobust.getValues(esol[i][j],E_ijl[i][j]);
        }
    }
    
    
    vector<vector<int>>vChains;
    vChains = GetChainsFrom1stStageSol(AdjacencyList,esol, Pairs, ChainLength);
    
    for (int i = 0; i < vChains.size(); i++){
        SolFirstStage.push_back(IndexGrandSubSol(vChains[i], vChains[i].size() - 1));
    }
    tTotal1stS += (clock() - tStart1stS)/double(CLOCKS_PER_SEC);

    
    //Call 2nd. stage
    tStart2ndS = clock();
    if (THP_Method != "PICEFBenders"){
        Chains2ndStage = Get2ndStageChains (SolFirstStage, RecoursePolicy);
    }
    Cycles2ndStage = Get2ndStageCycles (SolFirstStage, RecoursePolicy);
    tTotalFindingCyCh+= (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
    if (THP_Method == "Covering" || THP_Method == "DoubleCovering"){
        GrandSubProbMaster(Cycles2ndStage,Chains2ndStage,SolFirstStage);
    }
    else if (THP_Method == "Benders"){
        GrandSubProbMaster2(Cycles2ndStage,Chains2ndStage,SolFirstStage);
    }
    else if (THP_Method == "BendersPICEF"){
        GrandSubProbMasterPICEF(Cycles2ndStage, Chains2ndStage, SolFirstStage);
    }
    else{
        cout << "Invalid Method";
    }
    
}
