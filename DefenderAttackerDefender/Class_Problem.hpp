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
#include <algorithm>
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

struct vChain{
    int vertex;
    bool FirstPass = true;
    vector<int>veci;
    vector<int>::iterator it = veci.begin();
    vector<int>::iterator itEnd = veci.end();
    vChain(int _vertex, vector<int>sol){vertex = _vertex, veci = sol;}
};
struct Chain{
    vector<vChain> Vnodes;
    double AccumWeight = 0;
    Chain(vChain v){Vnodes.push_back(v);}
};
struct KEPSol{
    vector<int> cycles;
    vector<int> chains;
    KEPSol(vector<int> v_cy, vector<int> v_ch){cycles = v_cy, chains = v_ch;}
    KEPSol(){;}
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
    IloInt ChainLength;//L
    IloInt Ite2ndS = 0;
    IloInt Iteration = 0;
    IloRangeArray ActiveGrandSubSol;
    IloRangeArray BoundObjective;
    IloNumArray2 WeightMatrix;
    IloNumArray2 AdjacencyList;//Successors
    IloNumArray PRAList;//Successors
    string FileName;
    fstream file;
    clock_t tStart2ndS;
    double tEnd2ndS;
    string status;
    string FolderName;
    string RecoursePolicy;
    string THP_Method;
    string WhereItisRun;
    map<pair<int,int>,double>Weights;
    map<int,vector<int>>CycleNode;
    vector<vector<int>>PredList;
    map<int, vector<pair<int,int>>>PredMap;
    vector<Cycles> ListCycles;
    
    
    //Functions
    Problem(string _FolderName, string _FileName, IloInt _cycleLength, IloInt _chainLength, string _RecoursePolicy, string _THP_Method, IloInt _VertexBudget, IloInt _ArcBudget, string _WhereItisRun);
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
    IloNum FPMIP_Obj;
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
    map<int,vector<int>>CycleNodeSPH;
    map<int,vector<int>>ChainNodeSPH;
    map<int,vector<int>>CycleNodeTPH;
    map<int,vector<int>>ChainNodeTPH;
    
    
    //Constraint and Column Generation Grand Problem
    vector<int> GrandSubOptSol; //cycles and chains
    vector<int> GrandSubOptScenario; //Interdected cycles and chains
    void AddNewColsConsGP();
    
    //Grand SubProblem
    IloModel GrandSubProb;
    IloCplex cplexGrandSubP;
    IloNumVar Beta;
    IloNumVar Excess;
    IloObjective ObjGrandSubP;
    IloArray<IloNumColumn> NewCycleTPH;
    IloInt MaxArcFailures;
    IloInt MaxVertexFailures;
    IloInt RepSolCounter = 1;
    IloNum RobustObjTPH = 0;
    IloNum SPMIP_Obj = 0;
    IloNumVarArray cyvar;
    IloNumVarArray chvar;
    IloNumVarArray2 arc;
    IloNumVarArray vertex;
    IloNumVarArray selvertex;
    IloRangeArray TheOneCC;
    IloRangeArray MakeOneFailGrandSubP;
    IloRangeArray BoundObjGrandSubP;
    IloRangeArray ActiveCCSubP_CY;
    IloRangeArray ActiveCCSubP_CH;
    IloRangeArray IfNoFailures_CY;
    IloRangeArray IfNoFailures_CH;
    IloRangeArray Once;
    IloRangeArray SelVert2ndPH;
    IloRangeArray ConsBeta;
    IloRangeArray AtLeastOneFails;
    IloNumArray vertex_sol;
    IloNumArray2 arc_sol;
    map<pair<int,int>,int> mapArcs;
    vector<Cycles> RepairedListCCs;
    vector<Cycles> RobustSolTHP;
    vector<int> GrandProbSol;
    vector<vChain> VertexinSolChain;
    map<pair<int,int>, vector<int>> ArcsinCyclesTHP;
    map<pair<int,int>, vector<int>> ArcsinChainsTHP;
    map<pair<int,int>, bool> FailedArcs;
    map<int, bool> FailedVertices;
    map<pair<int,int>, bool> OptFailedArcs;
    map<int, bool> OptFailedVertices;
    map<int,int>Cycles2ndTo3rd;
    map<int,int>Chains2ndTo3rd;
    map<int,int>Cycles3rdTo2nd;
    map<int,int>Chains3rdTo2nd;
    map<int,vector<int>> CycleNodeGSP;
    
    void GrandSubProbMaster(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<IndexGrandSubSol>&SolFirstStage);
    void GrandSubProbRoutine();
    vector<Cycles> BackRecoursePolicy(vector<int>&vinFirstStageSol);
    vector<Cycles> AmongPolicy(vector<int>&vinFirstStageSol);
    vector<Cycles> AllPolicy(vector<int>&vinFirstStageSol);
    void AddNewColsConsGSP(vector<Cycles>& RepairedSol);
    void InitializeVertexinSolChain(vector<int>&ListVertices,vector<vChain>& VertexinSolChain);
    vector<Chain>FindChains(vector<vChain>& VertexinSolChain, vector<int>& vinFirstStageSol);
    vector<Chain> Get2ndStageChains (vector<IndexGrandSubSol>& GrandProbSol, string policy);
    vector<Cycles> Get2ndStageCycles (vector<IndexGrandSubSol>& GrandProbSol, string policy);
    void SampleCols2ndStage(vector<Chain>& Chains, vector<Cycles>&Cycles, vector<IndexGrandSubSol>&SolFirstStage);
    vector<int>GetSelVertices(vector<IndexGrandSubSol>&SolFirstStage);
    
    
    vector<int> Complete_ActiveCCSubP_LB(vector<int>PosNewCycles);
    void UpdateSNPSol(IloNumArray& r_sol, IloNum GrandSubObj);
    void PrintSolSNP(IloNumArray vertex_sol, IloNumArray2 arc_sol);
    void FillRobustSolTHP();
    
    //Third Phase
    IloObjective ObjTHP;
    vector<Cycles>Cycles2ndStage;
    vector<Chain>Chains2ndStage;
    void THPMIP(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<int>&ListSelVertices);
    IloNum VI_I = 0;
    //Model
    IloModel mTHPMIP;
    IloCplex cplexmTHPMIP;
    //Decision variables
    IloNumVarArray tcyvar;
    IloNumVarArray tchvar;
    IloNumVarArray Create_tcyvar(const char* prefix, vector<Cycles>&Cycles2ndStage);
    IloNumVarArray Create_tchvar(const char* prefix, vector<Chain>&Chains2ndStage);
    map<int,bool> GetUB_tcyvar(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices);
    map<int,bool> GetUB_tchvar(map<pair<int,int>, bool>&FailedArcs, map<int, bool>& FailedVertices);
    //Constraints
    IloRangeArray DisjointTHPArray;
    IloRangeArray DisjointTHP(IloNumVarArray& tcyvar, IloNumVarArray& tchvar);
    //Objective
    IloNum TPMIP_Obj = 0;
    IloNum LOWEST_TPMIP_Obj = INT_MAX;
    IloNumArray vals;
    int Cyclenewrow2ndPH = 0;
    int Chainnewrow2ndPH = 0;
    vector<int> tcysol3rd;
    vector<int> tchsol3rd;
    IloRange NewIloRangeCY3rd;
    IloRange NewIloRangeCH3rd;
    //Solution
    vector<Cycles>Cycles3rdSol;
    vector<Chain>Chains3rdSol;
    IloNumArray cyvar_sol2nd;
    IloNumArray chvar_sol2nd;
    //Algorithm
    void AddNewCols3rdTo2nd (IloNumArray tcysol, IloNumArray tchsol, map<int,int>& Cycles2ndTo3rd, map<int,int>& Chains2ndTo3rd, map<int,int>& Cycles3rdTo2nd, map<int,int>& Chains3rdTo2nd, int& Cyclenewrow2ndPH, int& Chainnewrow2ndPH, vector<int>& tcysol3rd, vector<int>& tchsol3rd, vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage);
    vector<int>ModifyOldActiveCCSubP_CY(int tOnecysol3rd, map<int,int>& Cycles2ndTo3rd, int& Cyclenewrow2ndPH, vector<Cycles>&Cycles2ndStage);
    vector<int>ModifyOldActiveCCSubP_CH(int tOnechsol3rd, map<int,int>& Chains2ndTo3rd, int& Chainnewrow2ndPH, vector<Chain>&Chains2ndStage);
    vector<int>ModifyOldSelectVex_CY(int tOnecysol3rd, vector<int>ListSelVertices, vector<Cycles>&Cycles2ndStage);
    vector<int>ModifyOldSelectVex_CH(int tOnechsol3rd, vector<int>ListSelVertices, vector<Chain>&Chains2ndStage);
    vector<int>ModifyOldActiveCCSubP_CYtoCH(int tOnecysol3rd, map<int,int>&Chains2ndTo3rd, int&Chainnewrow2ndPH, vector<Cycles>&Cycles2ndStage);
    vector<int>ModifyOldActiveCCSubP_CHtoCY(int tOnechsol3rd, map<int,int>& Cycles2ndTo3rd, int& Cyclenewrow2ndPH, vector<Chain>&Chains2ndStage);
    void selOnce(map<int,vector<int>>CycleNodeSPH, map<int,vector<int>>ChainNodeSPH);
    //void selOnceCH(map<int,vector<int>>CycleNodeSPH, map<int,vector<int>>ChainNodeSPH);
    void IfNoFailuresCY(map<int,int>& Cycles2ndTo3rd, int i);
    void IfNoFailuresCH(map<int,int>& Chains2ndTo3rd, int i);
    void GetNewIloRangeCY3rd(int tOnecysol3rd, vector<Cycles>&Cycles2ndStage);
    void GetNewIloRangeCH3rd(int tOnecysol3rd, vector<Chain>&Chains2ndStage);
    void GetNewBetaCut(IloNum TPMIP_Obj, map<pair<int,int>, bool> FailedArcs, map<int, bool> FailedVertices);
    void GetNoGoodCut(map<pair<int,int>, bool>& FailedArcs, map<int, bool>& FailedVertices);
    void GetAtLeastOneFails(vector<Cycles>&Cycles3rdSol, vector<Chain>&Chains3rdSol);
    void GetScenario(IloNumArray2& arc_sol, IloNumArray& vertex_sol);
    void Get3rdStageSol(vector<Cycles>&Cycles3rdSol, vector<Chain>&Chains3rdSol, IloNumArray& cyvar_sol3rd, IloNumArray& chvar_sol3rd);
    IloExpr GetObjTPH(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, string& TPH_Method);
    bool ThisWork(IloNumArray& tcysol, IloNumArray& tchsol, vector<int>&ListSelVertices, bool RemoveExcess);
        
    //Literature method
    void ROBUST_KEP();
    void GrandSubProbMaster2(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<IndexGrandSubSol>&SolFirstStage);
    bool Literature(IloNumArray& tcysol, IloNumArray& tchsol);
    void SampleCols2ndStage2(vector<Chain>& Chains, vector<Cycles>&Cycles, vector<IndexGrandSubSol>&SolFirstStage);
    vector<KEPSol>KEPSols2ndStage;
    IloRangeArray vBoundConstraint;
    IloRangeArray vFailedMatches;
    void Const11b(vector<KEPSol>&KEPSols2ndStage);
    void Const11c(vector<KEPSol>&KEPSols2ndStage);
    IloExpr exprVxtArcsCY;
    IloExpr exprVxtArcsCH;
    IloExpr exprBound;
    
    void Print2ndStage();
    
    void HeadingCF();
    void PrintCF();
private:
};


int checkIfCCisNew(vector<int>v, vector<IndexGrandSubSol>&Sol);
#endif /* Class_Problem_hpp */
