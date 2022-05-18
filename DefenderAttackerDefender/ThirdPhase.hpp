//
//  ThirdPhase.hpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-11-11.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef ThirdPhase_hpp
#define ThirdPhase_hpp

#include <stdio.h>
#include "Class_Problem.hpp"
#include "GrandSubproblem.hpp"
#include "GrandProblem.hpp"
#include "M-PICEF.hpp"


bool sortdouble(double& c1, double& c2);
bool sortCovElms(pair< pair<int,int>, double>& c1, pair< pair<int,int>, double>& c2);
bool sortduals(pair<int,double>& c1, pair<int,double>& c2);
bool sortPairs(pair<int,int>& c1, pair<int,int>& c2);
#endif /* ThirdPhase_hpp */
