//
//  Class_Problem.hpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-15.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef Class_Problem_hpp
#define Class_Problem_hpp

#include <stdio.h>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <map>

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2>  IloNumVarArray3;

class Cycles{
private:
    vector<int> _cycle;
    double _weight;
    int _HowMany;
public:
    Cycles(){}
    Cycles(vector<int> cycle, double weight){_cycle = cycle, _weight = weight;}
    double get_w(){return _weight;}
    vector<int> get_c(){return _cycle;}
    void set_Many(int HM){_HowMany = HM;}
    int get_Many(){return _HowMany;}
};

class IndexGrandSubSol{
private:
    vector<int> _GrandSubSol;
    double _weight;
    vector<int> _iteration;

public:
    IndexGrandSubSol(vector<int> GrandSubSol, double weight){_GrandSubSol = GrandSubSol, _weight = weight;}
    IndexGrandSubSol(vector<int> GrandSubSol, double weight, int ite){_GrandSubSol = GrandSubSol, _weight = weight, _iteration.push_back(ite);}
    void set_ite(int i){_iteration.push_back(i);}
    double get_ite(int i){return _iteration[i];}
    vector<int> get_itev(){return _iteration;}
    double get_w(){return _weight;}
    vector<int> get_cc(){return _GrandSubSol;}
};

class Problem{
public:
    IloEnv env;
    IloModel GrandProb;
    IloModel mPICEF;
    IloCplex cplexGrandP;
    IloNumColumn cgrandpsol;
    IloNumColumn cgrandpsolsce;
    IloNumColumn cgrandpparamsce;
    IloInt Pairs;
    IloInt NDDs;//altruists
    IloInt Nodes;//altruists + pairs
    IloInt NumArcs;
    IloInt CycleLength;//K
    IloInt ChainLength;//F
    IloInt Iteration = 0;
    IloRangeArray ActiveGrandSubSol;
    IloRangeArray BoundObjective;
    IloNumArray2 WeightMatrix;
    IloNumArray2 AdjacencyList;//Successors
    string FileName;
    string FolderName;
    map<pair<int,int>,double>Weights;
    map<int,vector<int>>CycleNode;
    vector<vector<int>>PredList;
    map<int, vector<pair<int,int>>>PredMap;
    vector<Cycles> ListCycles;
    int Hola;
    
    
    //Functions
    Problem(string _FolderName, string _FileName, IloInt _cycleLength, IloInt _chainLength);
    int Reading();
   
    //M-PICEF
    void M_PICEF();
    IloCplex cplexmPICEF;
    IloNumVarArray3 cyarc;
    IloNumVarArray3 charc;
    IloNumVarArray3 vrole;
    IloNumVarArray3 arole;
    IloNumVarArray2 arcori;
    IloNumVarArray2 aux;
    vector<int>fvs;
    
    void ModifyAdjacencyList(vector<int>v);
    IloNumVarArray3 CreateVar3Idx(IloInt maxindx, const char* prefix, vector<int> dist);// AdjacencyList, K or L
    IloNumVarArray3 CreateVar3Role(vector<int>fvs, IloNumArray2 g);
    IloNumVarArray2 CreateVarRole(vector<int>fvs); //Feedback Vertex Set
    IloNumVarArray3 CreateVar3ARole(vector<int>fvs, IloNumArray2 g);
    IloNumVarArray2 CreateVarArc2(const char* prefix, int upperbound); //AdjacencyList
    
    IloRangeArray Const1a();
    IloRangeArray Const1b();
    IloRangeArray Const1c();
    IloRangeArray Const1d();
    IloRangeArray Const1e();
    IloRangeArray Const1f();
    IloRangeArray Const1g();
    IloRangeArray Const1h();
    IloRangeArray Const1i();
    IloRangeArray Const1j();
    IloRangeArray Const1k();
    IloRangeArray Const1l();
    IloRangeArray Const1m();
    IloRangeArray Const1n();
    IloRangeArray ctest();
    
    //Get Objective
    IloExpr GetObjMPICEF();
    
    ////Dijkstra
    int minDistance(vector<int> dist, vector<bool> sptSet);
    void printSolution(vector<int> dist, int n);
    vector<int> dijkstra(IloNumArray2 graph, int src);
    void distCycles(IloNumArray2 graph);
    vector<int> distFor;
    vector<int> distBack;
    vector<vector<int>> distFVS;
    vector<vector<int>> distFVS_to_Pairs;
    vector<vector<int>> distPairs_to_FVS;
    map<pair<int,int>,vector<int>> ArcFvs;
    
    ////////Cycle Formulation/////
    IloNum Gap;
    IloNum SolTime;
    IloNum CycleSearchTime;
    IloNum NCycles;
    IloNum CFObj;
    IloNum BestObj;
    IloNumVarArray z;
    IloNumVar eta;
    void MainCycleFinder();
    void KEP_CycleFormulation();
    vector<Cycles> SubCycleFinder (IloEnv env, IloNumArray2 AdjaList, IloInt origin);
    void SetName1Index(IloNumVar& var, const char* prefix, IloInt i);
    IloBool IsxinStack (IloInt test, vector<int>& xinTrial);
    IloNum PathWeight (vector<int>& Stack);
    
    ////////Cycle Formulation Third Phase/////
    IloModel ThirdPH;
    IloCplex cplexThirdPH;
    vector<Cycles> CFThirdPhase(vector<Cycles>Input, map<int,vector<int>>&CycleMap);
    map<int,vector<int>>CycleNodeTPH;
    
    
    //Constraint and Column Generation Grand Problem
    vector<int> GrandSubOptSol; //cycles and chains
    vector<int> GrandSubOptScenario; //Interdected cycles and chains
    void AddNewColsConsGP();
    
    //Grand SubProblem
    IloModel GrandSubProb;
    IloCplex cplexGrandSubP;
    IloNumVar Beta;
    IloObjective ObjGrandSubP;
    IloArray<IloNumColumn> NewCycleTPH;
    IloInt MaxArcFailures = 2;
    IloInt MaxVertexFailures = 1;
    IloInt RepSolCounter = 1;
    IloNum RobustObjTPH = 0;
    IloNumVarArray r;
    IloNumVarArray2 arc;
    IloNumVarArray vertex;
    IloRangeArray TheOneCC;
    IloRangeArray MakeOneFailGrandSubP;
    IloRangeArray BoundObjGrandSubP;
    IloRangeArray ActiveCCSubP_LB;
    IloRangeArray ActiveCCSubP_UB;
    vector<Cycles> RepairedListCCs;
    vector<Cycles> RobustSolTHP;
    vector<IndexGrandSubSol> GrandSubSolSet;
    map<int,vector<int>> CycleNodeGSP;
    void GrandSubProbMaster();
    void GrandSubProbRoutine();
    vector<Cycles> BackRecoursePolicy(IloNumArray& vertex_sol);
    vector<Cycles> AmongPolicy(IloNumArray& vertex_sol);
    vector<Cycles> AllPolicy(IloNumArray& vertex_sol);
    void AddNewColsConsGSP(vector<Cycles>& RepairedSol);
    
    
    vector<int> Complete_ActiveCCSubP_LB(vector<int>PosNewCycles);
    void UpdateSNPSol(IloNumArray& r_sol, IloNum GrandSubObj);
    void PrintSolSNP(IloNumArray vertex_sol, IloNumArray2 arc_sol);
    void FillRobustSolTHP();

    void HeadingCF();
    void PrintCF();
private:
};


int checkIfCCisNew(vector<int>v, vector<IndexGrandSubSol>&Sol);
#endif /* Class_Problem_hpp */
