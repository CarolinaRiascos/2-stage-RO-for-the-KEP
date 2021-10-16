//
//  Dijkstra.cpp
//  DefenderAttackerDefender
//
//  Created by Lizeth Carolina Riascos Alvarez on 2021-10-06.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "Dijkstra.hpp"
// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
bool isAdjacent(IloNumArray row, int v);
bool fsptSet(vector<bool> vector, int v);
int Problem::minDistance(vector<int> dist, vector<bool> sptSet){
    // Initialize min value
    int min = INT_MAX, min_index;
 
    for (int v = 0; v < AdjacencyList.getSize(); v++)
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
 
    return min_index;
}
 
// A utility function to print the constructed distance array
void Problem::printSolution(vector<int> dist, int n)
{
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < AdjacencyList.getSize(); i++)
        printf("%d \t\t %d\n", i, dist[i]);
}
 
// Function that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
vector<int> Problem::dijkstra(IloNumArray2 graph, int src)
{
    vector<int> dist(graph.getSize()); // The output array.  dist[i] will hold the shortest
    // distance from src to i
 
    vector<bool> sptSet(graph.getSize()); // sptSet[i] will be true if vertex i is included in shortest
    // path tree or shortest distance from src to i is finalized
 
    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < graph.getSize(); i++)
        dist[i] = 2*Pairs, sptSet[i] = false;
 
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    int v; bool infsptSet = false;
    // Find shortest path for all vertices
    for (int count = 0; count < graph.getSize(); count++) {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in the first iteration.
        int u = minDistance(dist, sptSet);
 
        // Mark the picked vertex as processed
        sptSet[u] = true;
 
        if (dist[u] != 2*Pairs){
            // Update dist value of the adjacent vertices of the picked vertex.
            for (v = 0; v < Pairs; v++){
                // Update dist[v] only if is not in sptSet, there is an edge from
                // u to v, and total weight of path from src to  v through u is
                // smaller than current value of dist[v]
                //infsptSet = fsptSet(sptSet, v);
                bool visAdja = isAdjacent(graph[u], v);
                if (sptSet[v] == false &&  visAdja == true && dist[u] + 1 < dist[v]){
                    dist[v] = dist[u] + 1;
                }
            }
        }
    }
    // print the constructed distance array
    //printSolution(dist, int(graph.getSize()));
    
    return dist;
}
bool isAdjacent(IloNumArray row, int v){
    for (int i = 0; i < row.getSize(); i++){
        if (row[i] == v + 1) return true;
    }
    return false;
}
bool fsptSet(vector<bool> vector, int v){
    if (vector[v] == v) return 1;

    return 0;
}
void Problem::distCycles(IloNumArray2 graph){
    vector<int>dist(AdjacencyList.getSize(), INT_MAX);
    distFor = vector<int>(Pairs, 2*Pairs);
    distBack = vector<int>(Pairs, 2*Pairs);
    distFVS = vector<vector<int>> (fvs.size());
    distFVS_to_Pairs = vector<vector<int>> (fvs.size());
    distPairs_to_FVS = vector<vector<int>> (fvs.size());
    
    for (int i = 0; i < distFVS.size(); i++){
        distFVS[i] = vector<int>(fvs.size(),0);
        distFVS_to_Pairs[i] = vector<int>(Pairs,0);
        distPairs_to_FVS[i] = vector<int>(Pairs,0);
    }
    
    //One-way: From a feedback vertex to a pair
    for (int i = 0; i < fvs.size(); i++){
        dist = dijkstra(AdjacencyList, fvs[i]);
        for(int j = 0; j < dist.size(); j++){
            if (fvs[i] != j){
                distFVS_to_Pairs[i][j] = dist[j];
                if (dist[j] < distFor[j]) distFor[j] = dist[j];
                int pos = posinFVS(fvs, j);
                if (pos != - 1){
                    distFVS[i][pos] += dist[j];
                }
            }
        }
    }
    
    //Return: From pair to feedback vertices
    for(int i = 0; i < Pairs; i++){
        int posGo = posinFVS(fvs, i);
        dist = dijkstra(AdjacencyList, i);
        for(int j = 0; j < dist.size(); j++){
            int pos = posinFVS(fvs, j);
            if (pos != -1 && i != j){
                distPairs_to_FVS[pos][i] = dist[j];
                if (dist[j] < distBack[i]) distBack[i] = dist[j];
                if (posGo != -1){
                    distFVS[pos][posGo] = dist[j];
                }
            }
        }
    }
    
//    for (int i = 0; i < Pairs; i++){
//        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
//            for (int h = 0; h < fvs.size(); h++){
//                if (distFVS_to_Pairs[h][i] + distPairs_to_FVS[h][AdjacencyList[i][j] - 1] + 1 <= CycleLength){
//                    ArcFvs[make_pair(i,j)].push_back(h);
//                }
//            }
//        }
//    }
    
    cout << endl;
}
