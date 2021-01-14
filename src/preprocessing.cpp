/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

#include<iostream>
#include<fstream>
#include<vector>
#include<tuple>
#include<queue>
#include<algorithm>
#include<cstring>
#include"debug.h"
#include"stp_reader.cpp"
#include"graph.h"
using namespace std;

void dijkstra(int src, int n, vector<vector<pair<int, int>>> adj, int dist[]){
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
				if (mat[u][v] == -1){
					cout<<"Erro : "<<v<<' '<<u<<endl;
					return -1;
				}
				cnt++;
			}
	return cnt;			
}

void preprocessing(Graph &G){
	int n = G.V;
	int m = G. E;
	vector< vector< pair<int, int> > > adj = G.adj;

	bool term[n + 1];
	for (int v = 1; v <= n; v++)
		term[v] = G.p[v] > 0;

	int dist[n + 1][n + 1];

	bool ok = 1;
	while (ok){
		// least cost test;
		do{
			do{
				ok = 0;
			
				for (int v = 1; v <= n; v++)
					dijkstra(v, n, adj, dist[v]);

				int lst, vtx, dst;
			
				for (int v = 1; v <= n; v++){
					lst = adj[v].size() - 1;
	
					for (int i = 0; i <= lst; i++){
						tie(vtx, dst) = adj[v][i];
				//		cout<<v<<' '<<vtx<<' '<<dst<<endl;
						if (dist[v][vtx] < dst){
							swap(adj[v][i], adj[v][lst]);
							lst--;
							if (vtx >= v)
								m--;
							ok = 1;
						}
					}
					while (adj[v].size() > lst + 1)
						adj[v].pop_back();	
				}
			} while(ok);
		
		// Degree 1 test	

			ok = 0;
			for (int v = 1; v <= n; v++)
				if (adj[v].size() == 1 and G.p[v] > adj[v][0].second){
					int u = adj[v][0].first;
				
					for (int i = 0; i < adj[u].size(); i++)
						if (adj[u][i].first == v){ // maybe inneficient, but ok
							swap(adj[u][i], adj[u][adj[u].size() - 1]);
							adj[u].pop_back();
							ok = 1;
							break;
						}

					adj[v].clear();	
					m--;
					n--;	
				}
		} while (ok);	


	// for (auto i : adj[16])
	// 	dbg(i.first);

		// Now, we are going to add some stronger reductions that 	
		// need tranformations
		for (int v = 1; v <= n; v++)
			if (!term[v] and adj[v].size() == 2){
				ok = 1;
				int a = adj[v][0].first;
				int b = adj[v][1].first;
				int cst = adj[v][0].second + adj[v][1].second;	
			
		// 		cout<<"found "<<v<<' '<<a<<' '<<b<<endl;
		// 			for (auto i : adj[16])
		// dbg(i.first);

				for (int i = 0; i < adj[a].size(); i++)
					if (adj[a][i].first == v){
						dbg(adj[a][i].first);
						dbg(adj[a][adj[a].size() - 1].first);
						swap(adj[a][i], adj[a][adj[a].size() - 1]);
						adj[a].pop_back();
		// 				chapa;chapa;
		// 					for (auto i : adj[16])
		// dbg(i.first);

						break;
					}
					chapa;
			
				for (int i = 0; i < adj[b].size(); i++)
					if (adj[b][i].first == v){
						swap(adj[b][i], adj[b][adj[b].size() - 1]);
						adj[b].pop_back();
						break;
					}
				adj[v].clear();
	

				for (int i = 0; i < adj[a].size(); i++){
					if (adj[a][i].first == b){
						adj[a][i].second = min(adj[a][i].second, cst);
						break;
					}
					else if (i == adj[a].size() - 1){
						adj[a].push_back({b, cst});						
						m++;
					}
				}
				for (int i = 0; i < adj[b].size(); i++){
					if (adj[b][i].first == a){
						adj[b][i].second = min(adj[b][i].second, cst);
						break;
					}
					else if (i == adj[b].size() - 1)
						adj[b].push_back({a, cst});						
				}
				n--;
				m -= 2;
			}
	}

	int sum = 0;
	for (int v = 1; v <= n; v++)
		sum += adj[v].size();
	dbg(sum);

	int ssum = 0;
	for (int v = 1; v <= n; v++)
		for (auto u : adj[v])
			ssum++;
	dbg(ssum);	


	dbg(check(G.V, adj));

	G.adj = adj;
	G.E = m;
	G.V = n;

}


int main(){
	string s = "../instances/PCSPG-JMP/P100.4.stp";
	//string s = "../instances/morete.stp";

	Graph G;
	dbg(STP_reader(s, G));
	dbg(G.E);
	dbg(G.V);

	preprocessing(G);
	
	dbg(G.E);
	dbg(G.V);
}