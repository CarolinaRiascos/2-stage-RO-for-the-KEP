//
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-16.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef GrandProblem_hpp
#define GrandProblem_hpp

//#include <stdio.h>
//#include <lemon/list_graph.h> //Clase para crear la lista, hay varias pero esta es rápida
//#include <lemon/dijkstra.h> //Clase para que usemos dijkstra
//#include <lemon/connectivity.h>
//#include <ilcplex/ilocplex.h>
#include "Class_Problem.hpp"

ILOSTLBEGIN
IloRangeArray CreateCons7b(IloEnv& env, IloNumVar& Z, NumVar2D& Y_ju, int Ulast, int Pairs, string name);
IloRangeArray CreateCons7c(IloEnv& env, NumVar2D& Y_ju, IloNumVarArray& X_c, NumVar3D& E_ijl, int Ulast, int Pairs, map<int,vector<int>> CycleNode, map<int,vector<pair<int,int>>> PredMap, string name);
IloRangeArray CreateCon7d(IloEnv& env, NumVar2D& Y_ju, NumVar2D& X_cu, NumVar4D& E_ijlu, int Ulast, int Pairs, map<int,vector<int>> CycleNode, int Lmax, map<int,vector<pair<int,int>>> PredMap, string name);
IloRangeArray CreateCon7g(IloEnv& env, NumVar2D& X_cu, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, int Pairs, map<int,vector<int>> CycleNode, map<int, vector<pair<int,int>>> PredMap, int Lmax, string name);
IloRangeArray CreateCon7h(IloEnv& env, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, int P, IloNumArray2 Adja, string name);
IloRangeArray CreateCon7i(IloEnv& env, NumVar3D& E_ijl, int Pairs, int Lmax, map<int, vector<pair<int,int>>> PredMap, IloNumArray2 Adja, string name);
IloRangeArray CreateCon7e(IloEnv& env, IloNumVarArray& var2D, NumVar3D& var3D, int Pairs, map<int,vector<int>> CycleN, map<int,vector<pair<int,int>>> PredMap, string name);
IloRangeArray CreateCon7f(IloEnv& env, NumVar3D& E_ijl, int P, IloNumArray2 Adja, string name);
IloRangeArray CreateCon7j(IloEnv& env, NumVar4D& E_ijlu, int Ulast, int Pairs, int Lmax, map<int, vector<pair<int,int>>> PredMap, IloNumArray2 Adja, string name);
void CreateCon7k(IloEnv& env, NumVar2D& X_cu, NumVar4D& E_ijlu, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, IloNumArray2 AdjaList, map<pair<int,int>,vector<int>> CycleArcs, map<int, vector<pair<int,int>>> PredMap, int Lmax, string name);
IloRangeArray CreateCon7l(IloEnv& env, NumVar2D& X_cu, NumVar4D& E_ijlu, IloNumVarArray& X_c, NumVar3D& E_ijl, vector<map<pair<int,int>, bool>>&scenarios, int Ulast, map<int,vector<int>>CycleNode, map<int, vector<pair<int,int>>>PredMap, int Pairs, int Lmax, string name);
vector<vector<int>> GetChainsFrom1stStageSol(IloNumArray2 AdjacencyList,IloNumArray3 ysol, int Pairs, int ChainLength);
void SetUB3DBin (IloEnv& env, NumVar3D& vars, vector<int>& distNDD, int P, int L);
void SetUB4DBin (IloEnv& env, NumVar4D& vars, vector<int>& distNDD, int P, int L, int u);
NumVar3D Create3DBin (IloEnv& env, IloNumArray2 Adja, int L, string name);
#endif /* KEP_Deterministic_hpp */
