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

typedef vector< vector<pair<int, long long>>> vvii;
typedef vector< vector<long long>> vvi;
typedef vector<long long> vi;

void dijkstra(int src, int n, vvii adj, vi &dist){
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
int check(int n, vvii adj){
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

#endif 