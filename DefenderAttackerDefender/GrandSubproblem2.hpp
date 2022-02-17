//
//  GrandSubproblem2.hpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-12-03.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef GrandSubproblem2_hpp
#define GrandSubproblem2_hpp

#include <stdio.h>
#include "Class_Problem.hpp"
#include "M-PICEF.hpp"
ILOSTLBEGIN

bool sortChains2(Chain& c1, Chain& c2);
bool sortCycles2(Cycles& c1, Cycles& c2);
bool CheckRepetition(vector<int>& v, map<int,bool>& map);
bool CheckRepetition2(vector<vChain>& v, map<int,bool>& map);
#endif /* GrandSubproblem2_hpp */
