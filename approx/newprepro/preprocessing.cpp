/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <chrono>

#include "stp_reader.h"
#include "graph.h"
// #include "pcdist.cpp"
#include "dsu.h"
#include "frac.h"

using namespace std;

typedef vector< vector<pair<int, frac>>> vvif;
typedef vector< vector<int>> vvi;
typedef vector< vector<frac>> vvf;

typedef vector<int> vi;

void dijkstra(int src, int n, vvif adj, vector<frac> &dist){
	for (int v = 1; v <= n; v++)
		dist[v] = frac(-1);
	
	priority_queue< pair<frac, int>, vector<pair<frac, int> >, greater<pair<frac, int>> > pq;
	pq.push({frac(0), src});

	int vtx;
	frac dst;
	while (!pq.empty()){
		tie(dst, vtx) = pq.top();
		pq.pop();
		
		frac mn1 = frac(-1);
		if (dist[vtx] != mn1)
			continue;
	
		dist[vtx] = dst;
		
		for (auto x : adj[vtx])
			if (dist[x.first] == frac(-1))
				pq.push({dst + x.second, x.first});
	} 
}

// Check presolved graph for incosistencies, if one is found it returns -1
// otherwise it returns 2|E|
int check(int n, vvif adj){
	vector<vector<int>> mat(n + 1, vector<int>(n + 1, -1));
	for (int v = 1; v <= n; v++)
		for (auto u : adj[v])
			mat[v][u.first] = 1;

	int cnt = 0;	
	for (int v = 1; v <= n; v++)
		for (int u = 1; u <= n; u++)
			if (mat[v][u] == 1){
				if (mat[u][v] == -1)
					return -1;
				cnt++;
			}
	return cnt;			
}

bool least_cost(int &n, int &m, vvif &adj){
	bool ok = 0;
	vvf dist(n + 1, vector<frac>(n + 1));

	for (int v = 1; v <= n; v++) // calcula distâncias (O(nmlog(n)))
		dijkstra(v, n, adj, dist[v]);

	int lst, vtx;
	frac dst;

	for (int v = 1; v <= n; v++){
		lst = adj[v].size() - 1;

		for (int i = 0; i <= lst; i++){
			tie(vtx, dst) = adj[v][i];
			if (dist[v][vtx] < dst){
				swap(adj[v][i], adj[v][lst]);
				lst--;
				adj[v].pop_back();
				if (vtx >= v) // conta só uma vez
					m--;
				ok = 1; // retirei algo
			}
		}
	}

	return ok;
}

bool degree_1(int &n, int &m, vvif &adj, Graph &G){
	bool ok = 0;
	for (int v = 1; v <= n; v++){
		frac aux = frac(G.p[v]);
		if (adj[v].size() == 1 and aux < adj[v][0].second){
			int u = adj[v][0].first;
		
			for (int i = 0; i < adj[u].size(); i++)
				if (adj[u][i].first == v){ // maybe inneficient, but ok
					swap(adj[u][i], adj[u][adj[u].size() - 1]);
					adj[u].pop_back();
					break;
				}

			adj[v].clear();	
			m--;
			ok = 1;
		}
	}		
	return ok;	
}


// Preprocessing function, if some error occour the return value is -1,
// othwewise is 0.
int preprocessing(Graph &G){
	int n = G.V;
	int m = G.E;

	vvif adj(n + 1, vector<pair<int, frac>>(n + 1));
	for (int v = 1; v <= n; v++)
		for (auto u : G.adj[v])
			adj[v].push_back({u.first, frac(u.second)});

	vector<bool> term(n + 1);
	for (int v = 1; v <= n; v++)
		term[v] = G.p[v] > 0;

	bool ok = 1;
	while (ok){
		ok = 0;
		// ok |=	pcdist(n, m, p, term, adj);
		ok |= least_cost(n, m, adj);
		ok |= degree_1(n, m, adj, G);	
	}	
	
	int sum = 0;
	for (int v = 1; v <= G.V; v++)
		sum += adj[v].size();

	if (sum != check(G.V, adj))
		return -1;

	for (int v = 1; v <= n; v++){
		G.adj[v].clear();
		for (auto u : adj[v])
			G.adj[v].push_back({u.first, u.second.a});
	}

	G.E = m;
	return 0;
}

#ifndef PCST

int main(){
	// string s = "../instances/mst.stp";
	// string s = "../instances/PCSPG-JMP/P400.4.stp";
	string s = "../instances/E/E01-A.stp";

	Graph G;
	
	if (STP_reader(s, G) != 0){
		cout<<"Arquivo não existe"<<endl;
		return 0;
	}

	preprocessing(G);
}

#endif
