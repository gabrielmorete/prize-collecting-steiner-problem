/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

#include<iostream>
#include<fstream>
#include<vector>
// #include"stp_reader.cpp"
// #include"graph.h"
// #include"debug.h"
using namespace std;

// Returns the number of edges removed
int least_cost_test(Graph &G){
	int inf, adj[G.V][G.V];

	inf = 0;
	for (int v = 1; v <= G.V; v++)
		for (auto u : G.adj[v])
			inf += u.second;
	
	for (int v = 0; v < G.V; v++)
		for (int u = 0; u < G.V; u++)
			adj[v][u] = inf;

	for (int v = 0; v < G.V; v++)
		adj[v][v] = 0;

	for (int v = 1; v <= G.V; v++)
		for (auto u : G.adj[v])
			adj[v - 1][u.first - 1] = u.second;
	
	for (int k = 0; k < G.V; k++) // floyd-warshall
		for (int v = 0; v < G.V; v++)
			for (int u = 0; u < G.V; u++)
				adj[v][u] = min(adj[v][u], adj[v][k] + adj[k][u]);

	vector<pair<int, int>> adj_aux[G.V + 1];
			
	int rmv = 0;		
	for (int v = 1; v <= G.V; v++)
		for (auto u : G.adj[v])
			if (u.second <= adj[v - 1][u.first - 1])
				adj_aux[v].push_back(u);  // Bad, fix later
			else
				rmv++;
	
	for (int v = 1; v <= G.V; v++)
		G.adj[v] = adj_aux[v];

	G.E -= rmv / 2;

	return rmv / 2;		
}


// int main(){
// 	//string s = "../instances/B/b03.stp";
// 	string s = "../instances/PCSPG-JMP/K100.2.stp";

// 	Graph G;
// 	dbg(STP_reader(s, G));
// 	dbg(G.E);
// 	dbg(least_cost_test(G));
// 	dbg(G.E);
// }