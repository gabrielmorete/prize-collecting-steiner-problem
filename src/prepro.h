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

int check(int n, vector< vector< pair<int, int> > > adj){ // checa consistencia da saida
	int mat[n + 1][n + 1];
	memset(mat, -1, sizeof mat);

	for (int v = 1; v <= n; v++)
		for (auto u : adj[v])
			mat[v][u.first] = 1;

	int cnt = 0;	
	for (int v = 1; v <= n; v++)
		for (int u = 1; u <= n; u++)
			if (mat[v][u] == 1){
				if (mat[u][v] == -1){
					cout<<"Erro : "<<v<<' '<<u<<endl;
					return -1;
				}
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

// fundir o vértice u ao vértice v ?????? check
int fuse(int v, int u, int n, vector<vector<pair<int, int>>> &adj){
	int dst[n + 1];
	
	memset(dst, -1, sizeof dst);

	for (auto x : adj[v])
		dst[x.first] = x.second;

	for (auto x : adj[u]){
		if (dst[x.first] == -1)
			dst[x.first] = x.second;
		else
			dst[x.first] = min(dst[x.first], x.second);
	}

	for (auto x : adj[v]){	
		int a = x.first;
		if (a != u)
			for (int i = 0; i < adj[a].size(); i++)
				if (adj[a][i].first == v){
					adj[a].erase(adj[a].begin() + i);
					continue;
				}
	}		
	for (auto x : adj[u]){	
		int a = x.first;
		if (a != v)
			for (int i = 0; i < adj[a].size(); i++)
				if (adj[a][i].first == u){
					adj[a].erase(adj[a].begin() + i);
					continue;
				}
	}

	int clr = adj[v].size() + adj[u].size() - 1;
	adj[v].clear();
	adj[u].clear();

	for (int a = 1; a <= n; a++)
		if (a != v and a != u and dst[a] != -1){
			adj[v].push_back({a, dst[a]});
			adj[a].push_back({v, dst[a]});
			clr--;
		}
	
	return clr;
}

#endif 