/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
   Auxiliary functions
*/

#ifndef PREPRO
#define PREPRO

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<tuple>
#include<queue>
#include<algorithm>
#include<cstring>
#include<cassert>
using namespace std;

typedef vector< vector<pair<int, int>>> vvii;
typedef vector< vector<int>> vvi;
typedef vector<int> vi;

void dijkstra(int src, int n, vector<vector<pair<int, int>>> adj, vector<int> &dist){
	for (int v = 1; v <= n; v++)
		dist[v] = -1;
	
	priority_queue< pair<int, int>, vector<pair<int, int> >, greater<pair<int, int>> > pq;
	pq.push({0, src});

	int vtx, dst;
	while (!pq.empty()){
		tie(dst, vtx) = pq.top();
		pq.pop();
		
		if (dist[vtx] != -1)
			continue;
		dist[vtx] = dst;
		
		for (auto x : adj[vtx])
			if (dist[x.first] == -1)
				pq.push({dst + x.second, x.first});
	} 
}

// Check presolved graph for incosistencies, if one is found it returns -1
// otherwise it returns 2|E|
int check(int n, vector< vector< pair<int, int> > > adj){
	int mat[n + 1][n + 1];
	memset(mat, -1, sizeof mat);

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

void erase_adj(int v, int u, vector< vector< pair<int, int> > > &adj){
	for (int i = 0; i < adj[u].size(); i++)
		if (adj[u][i].first == v){
			adj[u].erase(adj[u].begin() + i);
			return;
		}
	assert(0);
}

bool add_adj(int a, int b, int cst, vector< vector< pair<int, int> > > &adj){
	bool is_ngh = 0;
	int ccst;
	for (auto x : adj[a])
		if (x.first == b){
			is_ngh = 1;
			ccst = x.second;
		}

	if (is_ngh == 0){ // add new edge
		adj[a].push_back({b, cst});
		adj[b].push_back({a, cst});
		return 1;
	}

	if (ccst <= cst)
		return 0;

	for (auto &x : adj[a])
		if (x.first == b){
			x.second = cst;
			break;
		}

	for (auto &x : adj[b])
		if (x.first == a){
			x.second = cst;
			break;
		}
	return 0;	
}

int get_min(int vtx, vvii &adj){
	int mn = adj[vtx][0].second;
	for (auto x : adj[vtx])
		mn = min(mn, x.second);
	return mn;
}

#endif 