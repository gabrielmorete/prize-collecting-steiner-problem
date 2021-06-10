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
#include "frac.h"
using namespace std;

typedef vector< vector<pair<int, frac>>> vvif;
typedef vector< vector<int>> vvi;
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
		
		if (dist[vtx] != frac(-1))
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

#endif 