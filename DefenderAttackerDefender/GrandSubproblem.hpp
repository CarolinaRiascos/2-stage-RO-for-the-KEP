//
//  GrandSubproblem.hpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-26.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef GrandSubproblem_hpp
#define GrandSubproblem_hpp

#include "Class_Problem.hpp"
#include "ThirdPhase.hpp"
#include "M-PICEF.hpp"
bool v2AlreadyinChain(vector<vChain> v1, int v2);
bool v2inFirstStageSol(vector<int>sol, int v);
bool sortNodes(Chain& c1, Chain& c2);
void FindNewNeighbor(vector<Chain>& PPChains);
ILOSTLBEGIN
#endif /* GrandSubproblem_hpp */
