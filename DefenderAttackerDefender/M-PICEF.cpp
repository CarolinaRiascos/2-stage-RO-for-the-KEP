//
//  M-PICEF.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-09-22.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "M-PICEF.hpp"



void Problem::M_PICEF(){
    //Create model
        mPICEF = IloModel (env);
        cplexmPICEF = IloCplex(mPICEF);
        cplexmPICEF.setParam(IloCplex::Param::TimeLimit, 250);
        cplexmPICEF.setParam(IloCplex::Param::Threads, 1);
    //Modify FVS
    //fvs = {1, 3, 6, 7, 18, 26, 29, 31, 34, 35, 39, 42, 59, 4, 28};
    //fvs = {0, 1, 2};
    fvs = {1, 7, 10, 11, 13, 15, 17, 23, 27, 32, 35, 41, 42, 49, 51, 53, 58, 62, 65, 66, 68, 76, 86, 91, 98, 100, 104, 108, 109, 112, 114, 122};
    
    //Dijkstra chains;
    vector<int>dist;
    vector<int>distNDD(AdjacencyList.getSize(), INT_MAX);
    for (int i = Pairs; i < AdjacencyList.getSize(); i++)
        dist = dijkstra(AdjacencyList, i);
        for(int j = 0; j < dist.size(); j++)
            if (dist[j] < distNDD[j]) distNDD[j] = dist[j];
    
    distCycles(AdjacencyList);
    
    //Create variables
    cyarc = CreateVar3Idx(CycleLength, "x", distNDD);
    charc = CreateVar3Idx(ChainLength, "y", distNDD);
    arcori = CreateVarArc2("s", int(fvs.size()));
    aux = CreateVarArc2("a", int(fvs.size()));
    vrole = CreateVar3Role(fvs, AdjacencyList);
    arole = CreateVar3ARole(fvs, AdjacencyList);
    
    //Constraint (1a)
    IloRangeArray Disjoint(env);
    Disjoint = Const1a();
    mPICEF.add(Disjoint);
    //Constraint (1b)
    IloRangeArray OneOutFromNDD(env);
    OneOutFromNDD = Const1b();
    mPICEF.add(OneOutFromNDD);
    //Constraint (1c)
    IloRangeArray OneOutFromFVS(env);
    OneOutFromFVS = Const1c();
    mPICEF.add(OneOutFromFVS);
    //Constraint (1d)
    IloRangeArray ChainPos(env);
    ChainPos = Const1d();
    mPICEF.add(ChainPos);
    //Constraint (1e)
    IloRangeArray CyclePos(env);
    CyclePos = Const1e();
    mPICEF.add(CyclePos);
    //Constraint (1f)
    IloRangeArray ArcBackToFVS(env);
    ArcBackToFVS = Const1f();
    mPICEF.add(ArcBackToFVS);
    //Constraint (1g)
    IloRangeArray OneRoleOnly(env);
    OneRoleOnly = Const1g();
    mPICEF.add(OneRoleOnly);
    //Constraint (1h)
    IloRangeArray BoundRoleToFlow(env);
    BoundRoleToFlow = Const1h();
    mPICEF.add(BoundRoleToFlow);
    //Constraint (1i)
    IloRangeArray SetFlowRoot(env);
    SetFlowRoot = Const1i();
    mPICEF.add(SetFlowRoot);
    //Constraint (1j)
    IloRangeArray EnforceRoot(env);
    EnforceRoot = Const1j();
    mPICEF.add(EnforceRoot);
//    //Constraint (1k)
//    IloRangeArray UBArcStart(env);
//    UBArcStart = Const1k();
//    mPICEF.add(UBArcStart);
    //Constraint (1l)
    IloRangeArray DominoEffect(env);
    DominoEffect = Const1l();
    mPICEF.add(DominoEffect);
    //Constraint (1m)
    IloRangeArray AuxConst(env);
    AuxConst = Const1m();
    mPICEF.add(AuxConst);
    //Constraint (1n)
//    IloRangeArray ExcludingRoles(env);
//    ExcludingRoles = Const1n();
//    mPICEF.add(ExcludingRoles);
    
    //Constraint (ctest)
    IloRangeArray test(env);
    test = ctest();
    //mPICEF.add(test);
    
    //Set Objective
    IloExpr obj(env,0);
    obj = GetObjMPICEF();
    mPICEF.add(IloMaximize(env, obj));
    
    //Solve M-PICEF
    cplexmPICEF.exportModel("MPICEF.lp");
    //cplexmPICEF.solve();
    
    IloNumArray3 xsol (env, Pairs);
    IloNumArray3 ysol (env, AdjacencyList.getSize());
    IloNumArray2 arcorisol (env, Pairs);
    IloNumArray2 auxsol (env, Pairs);
    IloNumArray2 rolesol (env, Pairs);

   if (cplexmPICEF.getStatus() == IloAlgorithm::Infeasible) {
       env.out() << "No solution" << endl;
   }
   else {
       //Retrieve solution

//       int Obj_MPICEF = cplexmPICEF.getObjValue();
//       env.out() << "M-PICEF Objective: " << Obj_MPICEF << endl;
       
       
//       for (int i = 0; i < xsol.getSize(); i++){
//           xsol[i] = IloNumArray2 (env, AdjacencyList[i].getSize());
//           for (int j = 0; j < xsol[i].getSize(); j++){
//               xsol[i][j] = IloNumArray(env);
//               cplexmPICEF.getValues(xsol[i][j],cyarc[i][j]);
//               for (int k = 0; k < CycleLength; k++){
//                   if (xsol[i][j][k] > 0){
//                       cout << cyarc[i][j][k].getName() << endl;
//                   }
//               }
//           }
//       }
//
//       for (int i = 0; i < ysol.getSize(); i++){
//           ysol[i] = IloNumArray2 (env, AdjacencyList[i].getSize());
//           for (int j = 0; j < ysol[i].getSize(); j++){
//               ysol[i][j] = IloNumArray(env, ChainLength);
//               cplexmPICEF.getValues(ysol[i][j],charc[i][j]);
//               for (int k = 0; k < ysol[i][j].getSize(); k++){
//                   if (ysol[i][j][k] > 0){
//                       cout << charc[i][j][k].getName() << endl;
//                   }
//               }
//           }
//       }
//
//       cout << "Origin: " << endl;
//       for (int i = 0; i < arcorisol.getSize(); i++){
//           arcorisol[i] = IloNumArray(env);
//           cplexmPICEF.getValues(arcorisol[i],arcori[i]);
//           for (int j = 0; j < arcorisol[i].getSize(); j++){
//               if (arcorisol[i][j] > 0){
//                   cout << arcori[i][j].getName() << ": " << arcorisol[i][j] << endl;
//               }
//           }
//       }

       vector<vector<int>>cysol;
       cysol.push_back(vector<int>{11, 112, 115,70, 11});
       cysol.push_back(vector<int>{24, 95, 27, 24});
       cysol.push_back(vector<int>{66, 17, 66});
       cysol.push_back(vector<int>{67, 4, 92, 12, 67});
       cysol.push_back(vector<int>{69, 108, 50, 75, 69});
       cysol.push_back(vector<int>{105, 49, 76, 85, 105});
       vector<vector<int>>chsol;
       chsol.push_back(vector<int>{129, 39, 52, 26});
       chsol.push_back(vector<int>{130, 14, 25, 59});
       chsol.push_back(vector<int>{131, 30, 57, 93});
       chsol.push_back(vector<int>{132, 122, 63, 36});
       chsol.push_back(vector<int>{133, 2, 55, 87});
       chsol.push_back(vector<int>{134, 6, 28, 33});
       
       vector<IndexGrandSubSol>solcy; vector<IndexGrandSubSol>solch;
       for (int i = 0; i < 6; i++){
           solcy.push_back(IndexGrandSubSol(cysol[i], cysol[i].size() - 1));
           solch.push_back(IndexGrandSubSol(chsol[i], chsol[i].size()));
       }
       
       vector<Cycles>Cycles2ndStage;
       vector<Chain>Chains2ndStage;
       Chains2ndStage = Get2ndStageChains (solch, RecoursePolicy);
       Cycles2ndStage = Get2ndStageCycles (solcy, RecoursePolicy);
       GrandSubProbMaster(Cycles2ndStage,Chains2ndStage);
   }
   
}
IloExpr Problem::GetObjMPICEF(){
    IloExpr obj(env,0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            if (i < Pairs){
                obj+= IloSum(charc[i][j]) + IloSum(cyarc[i][j]);//
            }
            else{
                obj+= IloSum(charc[i][j]);
            }
        }
    }
    return obj;
}
void Problem::ModifyAdjacencyList(vector<int>fvs){
    
    IloNumArray2 NewAdjaList(env, Pairs + NDDs + fvs.size());
    for (int i = 0; i < NewAdjaList.getSize(); i++){
        NewAdjaList[i] = IloNumArray (env);
    }
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        int row = i;
        if (i >= Pairs){
            row = i + fvs.size();
        }
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            int idx = posinFVS(fvs, AdjacencyList[i][j] - 1);
            if (idx != -1){
                NewAdjaList[row].add(Pairs + idx + 1);
            }
            else{
                NewAdjaList[row].add(AdjacencyList[i][j]);
            }
        }
    }
    for (int i = 0; i < NewAdjaList.getSize(); i++){
        cout <<  endl << i + 1 << ": " ;
        for (int j = 0; j < NewAdjaList[i].getSize(); j++){
            //cout << NewAdjaList[i][j] << " " << "\t";
        }
    }
    AdjacencyList.clear();
    AdjacencyList = NewAdjaList;
    

}
IloRangeArray Problem::ctest(){
    IloRangeArray test(env);
    IloExpr expr (env,0);
    
    for (int j = 0; j < AdjacencyList[fvs[0]].getSize(); j++){
        expr+= cyarc[fvs[0]][j][0];
    }
    
    string name = "ctest";
    const char* cName = name.c_str();
    test.add(IloRange(env, 1, cyarc[0][0][0], IloInfinity, cName));
    
    return test;
}
IloRangeArray Problem::Const1a(){
    IloRangeArray Disjoint(env);
    
    for (int i = 0; i < Pairs; i++){
        IloExpr expr (env,0);
        for (int j = 0; j < PredMap[i].size(); j++){
            if (PredMap[i][j].first < Pairs){
                expr +=  IloSum(cyarc[PredMap[i][j].first][PredMap[i][j].second]) + IloSum(charc[PredMap[i][j].first][PredMap[i][j].second]);
            }
            else{
                expr+= IloSum(charc[PredMap[i][j].first][PredMap[i][j].second]);
            }
        }
        string name = "Ca." + to_string(i + 1);
        const char* cName = name.c_str();
        Disjoint.add(IloRange(env, -IloInfinity, expr, 1, cName));
    }
    return Disjoint;
}
IloRangeArray Problem::Const1b(){
    IloRangeArray OneOutFromNDD(env);
    
    for (int i = Pairs; i < AdjacencyList.getSize(); i++){
        IloExpr expr (env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            expr += charc[i][j][0];
        }
        string name = "Cb." + to_string(i + 1);
        const char* cName = name.c_str();
        OneOutFromNDD.add(IloRange(env, -IloInfinity, expr, 1, cName));
    }
    return OneOutFromNDD;
}
IloRangeArray Problem::Const1c(){
    IloRangeArray OneOutFromFVS(env);
    
    for (int i = 0; i < Pairs; i++){
        IloExpr expr (env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            expr += IloSum(cyarc[i][j]);
        }
        string name = "Cc."  + to_string(i + 1);
        const char* cName = name.c_str();
        OneOutFromFVS.add(IloRange(env, -IloInfinity, expr, 1, cName));
    }
    
    return OneOutFromFVS;
}
IloRangeArray Problem::Const1d(){
    IloRangeArray ChainPos(env);

    for (int i = 0; i < Pairs; i++){
        for (int k = 0; k < ChainLength - 1; k++){
            IloExpr exprIn (env,0);
            IloExpr exprOut (env,0);
            for (int j = 0; j < PredMap[i].size(); j++){
                if (PredMap[i][j].first < Pairs){
                    exprIn += charc[PredMap[i][j].first][PredMap[i][j].second][k];
                }
                else if (k == 0){
                    exprIn += charc[PredMap[i][j].first][PredMap[i][j].second][0];
                }
            }
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                exprOut += charc[i][j][k + 1];
            }
            //cout << exprIn << " > " << exprOut << endl;
            string name = "Cd."  + to_string(i + 1);
            const char* cName = name.c_str();
            ChainPos.add(IloRange(env, 0, exprIn - exprOut, IloInfinity, cName));
        }
    }
    return ChainPos;
}
IloRangeArray Problem::Const1e(){
    IloRangeArray CyclePos(env);
    
    for (int i = 0; i < Pairs; i++){
        for (int k = 0; k < CycleLength - 1; k++){
            IloExpr exprIn (env,0);
            IloExpr exprOut (env,0);
            for (int j = 0; j < PredMap[i].size(); j++){
                if (PredMap[i][j].first < Pairs){
                    exprIn += cyarc[PredMap[i][j].first][PredMap[i][j].second][k];
                }
            }
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                exprOut += cyarc[i][j][k + 1];
            }
            string name = "Ce."  + to_string(i + 1);
            const char* cName = name.c_str();
            //cout << exprIn << " >= " << exprOut << endl;
            CyclePos.add(IloRange(env, 0, exprIn - exprOut, IloInfinity, cName));
        }
    }
    return CyclePos;
    
}
//IloRangeArray Problem::Const1f(){
//    IloRangeArray ArcBackToFVS(env);
//
//    for (int i = Pairs; i < Pairs + fvs.size(); i++){
//        IloExpr exprIn (env,0);
//        for (int j = 0; j < PredMap[i].size(); j++){
//            if (PredMap[i][j].first < Pairs){
//                for (int k = 1; k < CycleLength; k++){
//                    exprIn += cyarc[PredMap[i][j].first][PredMap[i][j].second][k];
//                }
//            }
//        }
//        //cout << exprIn << " = " << vrole[i - Pairs].getName() << endl << endl;
//        string name = "Cf."  + to_string(i + 1);
//        const char* cName = name.c_str();
//        ArcBackToFVS.add(IloRange(env, 0, exprIn - vrole[i - Pairs], 0, cName));
//    }
//
//    return ArcBackToFVS;
//}
//IloRangeArray Problem::Const1g(){
//    IloRangeArray NoCompIfnoFVSin1(env);
//    for (int i = 0; i < fvs.size(); i++){
//        IloExpr expr (env, 0);
//        for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
//            expr+= cyarc[fvs[i]][j][0];
//        }
//        //cout << vrole[i].getName()  << " = " << expr << endl;
//        string name = "Cg."  + to_string(fvs[i] + 1);
//        const char* cName = name.c_str();
//        NoCompIfnoFVSin1.add(IloRange(env, 0, vrole[i] - expr, 0, cName));
//    }
//
//    return NoCompIfnoFVSin1;
//}

IloRangeArray Problem::Const1f(){
    IloRangeArray ArcBackToFVS(env);

    for (int i = 0; i < fvs.size(); i++){
        IloExpr exprIn (env,0);
        IloExpr exprOut (env,0);
        for (int j = 0; j < PredMap[fvs[i]].size(); j++){
            if (PredMap[fvs[i]][j].first < Pairs){
                for (int k = 1; k < CycleLength; k++){
                    exprIn += cyarc[PredMap[fvs[i]][j].first][PredMap[fvs[i]][j].second][k];
                }
            }
        }
        for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
            exprOut += cyarc[fvs[i]][j][0];
        }
        //cout << exprIn << " > " << exprOut << endl << endl;
        string name = "Cf."  + to_string(i + 1);
        const char* cName = name.c_str();
        ArcBackToFVS.add(IloRange(env, 0, exprIn - exprOut, IloInfinity, cName));
    }

    return ArcBackToFVS;
}
IloRangeArray Problem::Const1g(){
    IloRangeArray OneRole(env);
    for (int u = 0; u < fvs.size(); u++){
        for (int v = 0; v < AdjacencyList[fvs[u]].getSize(); v++){
            string name = "Cg."  + to_string(fvs[u] + 1) + "." + to_string(AdjacencyList[fvs[u]][v]);
            const char* cName = name.c_str();
            //cout << IloSum(vrole[u][v]) << endl;
            OneRole.add(IloRange(env, -IloInfinity, IloSum(vrole[u][v]) - IloSum(cyarc[fvs[u]][v]), 0, cName));
        }
    }

    return OneRole;
}
IloRangeArray Problem::Const1h(){
    IloRangeArray RoleisToStart(env);
    int u,i;
    for (u = 0; u < fvs.size(); u++){
        for (int v = 0; v < AdjacencyList[fvs[u]].getSize(); v++){
            for (i = 0; i < fvs.size(); i++){
                if (i >= u){
                    IloExpr expr (env, 0);
                    for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
                        if (i == u){
                            expr+= cyarc[fvs[i]][j][0];
                        }
                        else{
                            if (distPairs_to_FVS[u][AdjacencyList[fvs[i]][j] - 1] + distPairs_to_FVS[i][AdjacencyList[fvs[u]][v] - 1] + 2 <= CycleLength){
                                expr+= cyarc[fvs[i]][j][0];
                            }
                        }
                    }
                    string name = "Ch."  + to_string(fvs[u] + 1) + "." + to_string(AdjacencyList[fvs[u]][v]);
                    const char* cName = name.c_str();
                    //cout << vrole[u][v][i].getName() << " <= " << expr << endl;
                    RoleisToStart.add(IloRange(env, -IloInfinity, vrole[u][v][i] - expr, 0, cName));
                }
            }
        }
    }
    return RoleisToStart;
}
bool checkIfArcFeasible(vector<vector<int>> distFVS_to_Pairs, vector<vector<int>> distPairs_to_FVS, pair<int,int> pair){
    
    
    
    return false;
}
IloRangeArray Problem::Const1i(){
    IloRangeArray SetFlowRoot(env);

    for (int u = 0; u < fvs.size(); u++){
        for (int v = 0; v < AdjacencyList[fvs[u]].getSize(); v++){
            IloExpr expr (env, 0);
            for (int i = 0; i < fvs.size(); i++){
                if (distFVS_to_Pairs[i][fvs[u]] + distPairs_to_FVS[i][AdjacencyList[fvs[u]][v] - 1] + 1 <= CycleLength){
                    expr += (i + 1)*vrole[u][v][i];
                }
            }
            string name = "Ci."  + to_string(fvs[u] + 1) + to_string(AdjacencyList[fvs[u]][v]);
            const char* cName = name.c_str();
            //cout << arcori[fvs[u]][v].getName() << " == " << expr << endl;
            SetFlowRoot.add(IloRange(env, 0, arcori[fvs[u]][v] - expr, 0, cName));
        }
    }
    
    for (int u = 0; u < Pairs; u++){
        for (int v = 0; v < AdjacencyList[u].getSize(); v++){
            if (inFVS(fvs, u) == false){
                IloExpr expr (env, 0);
                for (int i = 0; i < fvs.size(); i++){
                    if (distFVS_to_Pairs[i][u] + distPairs_to_FVS[i][AdjacencyList[u][v] - 1] + 1 <= CycleLength){
                        expr += (i + 1)*arole[u][v][i];
                    }
                }
                string name = "Ci."  + to_string(u + 1) + to_string(AdjacencyList[u][v]);
                const char* cName = name.c_str();
                //cout << arcori[fvs[u]][v].getName() << " == " << expr << endl;
                SetFlowRoot.add(IloRange(env, 0, arcori[u][v] - expr, 0, cName));
            }
        }
    }
    
    return SetFlowRoot;
}
//IloRangeArray Problem::Const1j(){
//    IloRangeArray SetFVSArcStart(env);
//
//    for (int i = 0; i < fvs.size(); i++){
//        IloExpr expr(env, 0);
//        int pos = posinFVS(fvs, fvs[i]);
////        for (int h = 0; h < fvs.size(); h++){
////            if (h != pos)  expr+= (h+1)*vrole[h][pos];
////        }
//        for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
//            string name = "Cj."  + to_string(fvs[i] + 1) + "." + to_string(AdjacencyList[fvs[i]][j]);
//            const char* cName = name.c_str();
//            //cout << arcori[fvs[i]][j].getName() << " = " << expr << "+" << fvs[i] + 1 << cyarc[fvs[i]][j][0].getName() << endl;
//            SetFVSArcStart.add(IloRange(env, -IloInfinity, arcori[fvs[i]][j] - (fvs.size() - i - 1)*( 1 - cyarc[fvs[i]][j][0]) - (i + 1), 0, cName));
//            SetFVSArcStart.add(IloRange(env, 0, arcori[fvs[i]][j] - (i + 1)*cyarc[fvs[i]][j][0], IloInfinity, cName));
//        }
//    }
//
//    return SetFVSArcStart;
//}
IloRangeArray Problem::Const1j(){
    IloRangeArray EnforceRoot(env);
    
    for (int i = 0; i < fvs.size(); i++){
        for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
            string name = "Cj."  + to_string(fvs[i] + 1) + to_string(AdjacencyList[fvs[i]][j]);
            const char* cName = name.c_str();
            //cout << cyarc[fvs[i]][j][0].getName() << " == " << vrole[i][j][i].getName() << endl;
            EnforceRoot.add(IloRange(env, 0, cyarc[fvs[i]][j][0] - vrole[i][j][i], 0, cName));
        }
    }
    
    return EnforceRoot;
}
IloRangeArray Problem::Const1k(){
    IloRangeArray UBArcStart(env);

    for (int i = 0; i < Pairs; i++){
        for (int j = 0; j  < AdjacencyList[i].getSize(); j++){
            string name = "Ck."  + to_string(i + 1) + "." + to_string(AdjacencyList[i][j]);
            const char* cName = name.c_str();
            //cout << fvs.size()*IloSum(cyarc[i][j]) << endl;
            UBArcStart.add(IloRange(env, -IloInfinity, arcori[i][j] - fvs.size()*IloSum(cyarc[i][j]), 0, cName));
        }
    }

    return UBArcStart;
}
//IloRangeArray Problem::Const1l(){
//    IloRangeArray DominoEffect(env);
//
//    for (int i = 0; i < Pairs; i++){
//        IloExpr expr(env, 0);
//        for (int j = 0; j < PredMap[i].size(); j++){
//            if (PredMap[i][j].first < Pairs){
//                expr+= arcori[PredMap[i][j].first][PredMap[i][j].second];
//            }
//        }
//        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
//            string name = "Cl."  + to_string(i + 1) + "." + to_string(AdjacencyList[i][j]);
//            const char* cName = name.c_str();
//            DominoEffect.add(IloRange(env, -IloInfinity, arcori[i][j] - expr , 0, cName));//+ aux[i][j]
//        }
//    }
//
//    return DominoEffect;
//}
IloRangeArray Problem::Const1l(){
    IloRangeArray DominoEffect(env);
    
    for(int i = 0; i < Pairs; i++){
        IloExpr exprIn (env,0);
        IloExpr exprOut (env,0);
        for (int j = 0; j < PredMap[i].size(); j++){
            if (PredMap[i][j].first < Pairs){
                exprIn += arcori[PredMap[i][j].first][PredMap[i][j].second];
            }
        }
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            exprOut += arcori[i][j];
        }
        string name = "Cl."  + to_string(i + 1);
        const char* cName = name.c_str();
        //cout << exprIn << " = " << exprOut << endl;
        DominoEffect.add(IloRange(env, 0, exprIn - exprOut, 0, cName));
    }
    
    return DominoEffect;
}
IloRangeArray Problem::Const1m(){
    IloRangeArray AuxConst(env);

    for (int i = 0; i < Pairs; i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            string name = "Cm."  + to_string(i + 1) + "." + to_string(AdjacencyList[i][j]);
            const char* cName = name.c_str();
            //cout << arcori[i][j] - fvs.size()*IloSum(cyarc[i][j]) << endl;
            AuxConst.add(IloRange(env, -IloInfinity, arcori[i][j] - fvs.size()*IloSum(cyarc[i][j]), 0, cName));
            //AuxConst.add(IloRange(env, -IloInfinity, aux[i][j] - arcori[i][j], 0));
        }
    }
    
    return AuxConst;
}
IloRangeArray Problem::Const1n(){
    IloRangeArray ExcludingRoles(env);

    for (int i = 0; i < fvs.size(); i++){
        for (int h = 0; h < fvs.size(); h++){
            IloExpr expr1(env,0);
            IloExpr expr2(env,0);
            if (i != h){
                for (int j = 0; j < AdjacencyList[fvs[i]].getSize(); j++){
                    expr1+= vrole[i][j][h];
                }
                for (int j = 0; j < AdjacencyList[fvs[h]].getSize(); j++){
                    expr2+= vrole[h][j][i];
                }
                //cout << expr1 + expr2 << endl << endl;
                string name = "Cn."  + to_string(fvs[i] + 1) + to_string(fvs[h] + 1);
                const char* cName = name.c_str();
                ExcludingRoles.add(IloRange(env, -IloInfinity, expr1 + expr2, 1, cName));
            }
        }
    }

    return ExcludingRoles;
}

IloNumVarArray3 Problem::CreateVar3Idx(IloInt maxindx, const char* prefix, vector<int> distNDD){
    int upto = 0;
    if (prefix == "x"){
        upto = Pairs;
    }else{upto = AdjacencyList.getSize();}
    IloNumVarArray3 var(env, upto);
    for (int i = 0; i < upto; i++){
        var[i] = IloNumVarArray2 (env, AdjacencyList[i].getSize());
        for (int j = 0; j < var[i].getSize(); j++){
            if (i < Pairs){
                var[i][j] = IloNumVarArray(env, maxindx, 0, 1, ILOINT);
            }
            else{
                var[i][j] = IloNumVarArray(env, 1, 0, 1, ILOINT);
            }
            for (int k = 0; k < var[i][j].getSize(); k++){
                SetName3(var[i][j][k], prefix, i + 1, AdjacencyList[i][j], k + 1);
                if (upto == AdjacencyList.getSize() && i < Pairs){
                    if (k < distNDD[i]) var[i][j][k].setUB(0);
                }
                if (upto == Pairs && k == 0 && inFVS(fvs,i) == false){
//                    if (distFor[i] + distBack[AdjacencyList[i][j] - 1] + 1 <= CycleLength){
//                        if (k < distFor[i]){
//                            if (inFVS(fvs,i) == false){
//                                var[i][j][k].setUB(0);
//                            }
//                            else{
//                                //cout << distFor[i] << endl;
//                                if (k != 0) var[i][j][k].setUB(0);
//                            }
//                        }
//                    }
//                    else{
//                        if (inFVS(fvs,i) == false){
//                            var[i][j][k].setUB(0);
//                        }
//                        else{
//                            if (k != 0) var[i][j][k].setUB(0);
//                        }
//                    }
                    var[i][j][k].setUB(0);
                }
                if (upto == AdjacencyList.getSize() && ChainLength <= 0) var[i][j][k].setUB(0);
                //cout << var[i][j][k].getName() << endl;
            }
        }
    }
    return var;
}
bool inFVS (vector<int> v, int vertex){
    for (int i = 0; i < v.size(); i++){
        if (vertex == v[i]) return true;
    }
    return false;
}
IloNumVarArray3 Problem::CreateVar3Role(vector<int>fvs, IloNumArray2 g){
    IloNumVarArray3 var(env, fvs.size());
    
    for (int i = 0; i < var.getSize(); i++){
        var[i] = IloNumVarArray2 (env, AdjacencyList[fvs[i]].getSize());
        for (int j = 0; j < var[i].getSize(); j++){
            var[i][j] = IloNumVarArray(env, fvs.size(), 0, 1, ILOINT);
            for (int k = 0; k < var[i][j].getSize(); k++){
                SetName3(var[i][j][k], "r", fvs[i] + 1, AdjacencyList[fvs[i]][j], fvs[k] + 1);
            }
        }
    }
    return var;
    
}
IloNumVarArray3 Problem::CreateVar3ARole(vector<int>fvs, IloNumArray2 g){
    IloNumVarArray3 var(env, AdjacencyList.getSize());
    
    for (int i = 0; i < var.getSize(); i++){
        if (inFVS(fvs, i) == false){
            var[i] = IloNumVarArray2 (env, AdjacencyList[i].getSize());
            for (int j = 0; j < var[i].getSize(); j++){
                var[i][j] = IloNumVarArray(env, fvs.size(), 0, 1, ILOINT);
                for (int k = 0; k < var[i][j].getSize(); k++){
                    SetName3(var[i][j][k], "r", i + 1, AdjacencyList[i][j], fvs[k] + 1);
                }
            }
        }
    }
    return var;
    
}
IloNumVarArray2 Problem::CreateVarRole(vector<int>fvs){
    IloNumVarArray2 var(env, fvs.size());
    for (int i = 0; i < var.getSize(); i++){
        var[i] = IloNumVarArray(env, fvs.size(), 0, 1, ILOINT);
        for (int j = 0; j < var[i].getSize(); j++){
            SetName2(var[i][j], "r", fvs[i] + 1, fvs[j] + 1);
            //cout << var[i][j] << endl;
        }
    }
    return var;
}
//IloNumVarArray Problem::CreateVarRole(vector<int>fvs){
//    IloNumVarArray var(env, fvs.size(), 0, 1, ILOINT);
//
//    for (int i = 0; i < var.getSize(); i++){
//        SetName(var[i], "r", fvs[i] + 1);
//    }
//
//    return var;
//}
IloNumVarArray2 Problem::CreateVarArc2(const char* prefix, int upperbound){
    IloNumVarArray2 var(env, Pairs);
    for (int i = 0; i < var.getSize(); i++){
        var[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, upperbound, ILOINT);
        for (int j = 0; j < var[i].getSize(); j++){
            SetName2(var[i][j], prefix, i + 1, AdjacencyList[i][j]);
        }
    }
    return var;
}
void SetName(IloNumVar& var, const char* prefix, IloInt i){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i);
    const char* varName = name.c_str();
    var.setName(varName);
}
void SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
void SetName3(IloNumVar& var, const char* prefix, IloInt e, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(e) + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
int posinFVS(vector<int>fvs, int lookfor){
    for (int i = 0; i < fvs.size(); i++){
        if (lookfor == fvs[i]) return i;
    }
    return -1;
}
