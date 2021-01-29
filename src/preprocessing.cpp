/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

<<<<<<< HEAD
#include"debug.h"
#include"prepro.h"
#include"stp_reader.h"
#include"graph.h"
#include"pcdist.cpp"
using namespace std;

bool least_cost(int &n, int &m, vvii &adj, vvi &dist){
	bool ok = 0;
	for (int v = 1; v <= n; v++) // calcula distâncias (O(nmlog(n)))
		dijkstra(v, n, adj, dist[v]);

	int lst, vtx, dst;

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

bool degree_1(int &n, int &m, vvii &adj, Graph &G){
	bool ok = 0;
	for (int v = 1; v <= n; v++)
		if (adj[v].size() == 1 and G.p[v] < adj[v][0].second){
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
	return ok;	
}

bool degree_2(int &n, int &m, int &en, vector<bool> term, vvii &adj){
	bool ok = 0;
	for (int v = 1; v <= n; v++)
		if (!term[v] and adj[v].size() == 2){
			ok = 1;

			int a = adj[v][0].first;
			int b = adj[v][1].first;
			int cst = adj[v][0].second + adj[v][1].second;

			erase_adj(v, a, adj);
			erase_adj(v, b, adj);
		
			adj[v].clear();
			en--;
			m -= 2 - add_adj(a, b, cst, adj);
	}
	return ok;
}

bool degree_3(int &n, int &m, int &en, vector<bool> &term, vvii &adj, vvi &dist){
	bool ok = 0;

	for (int v = 1; v <= n; v++)
		if (!term[v] and adj[v].size() == 3){

			int a = adj[v][0].first, va = adj[v][0].second;
			int b = adj[v][1].first, vb = adj[v][1].second;
			int c = adj[v][2].first, vc = adj[v][2].second;
			int cst = va + vb + vc;
			int dst = min({dist[b][a] + dist[a][c], dist[a][b] + dist[b][c], dist[a][c] + dist[c][b]});

			if (dist[a][b] == -1 or dist[a][c] == -1 or dist[b][c] == -1)
				assert(0);

			if (dst < cst){
				ok = 1;
				adj[v].clear();
				en--;
				m -= 3;
				
				erase_adj(v, a, adj);
				erase_adj(v, b, adj);
				erase_adj(v, c, adj);

				m += add_adj(a, b, va + vb, adj);
				m += add_adj(a, c, va + vc, adj);
				m += add_adj(b, c, vb + vc, adj);
			}
		}	
	return ok;	
}

bool min_adj(int &n, int &m, int &en, vector<bool> &term, vvii &adj, Graph &G){
	bool ok = 0;
	for (int v = 1; v <= n; v++){
		if (!term[v] or adj[v].empty())
			continue;
	
		int cmn = adj[v][0].second;
		for (auto x : adj[v])
			cmn = min(cmn, x.second);

		for (auto x : adj[v]){
			int u = x.first;
			int cst = x.second;

			if (!term[u] or cmn != cst)
				continue;

			if (min(G.p[v], G.p[u]) > cst){
				m -= fuse(v, u, n, adj);
				G.p[v] = G.p[v] + G.p[u] - cst;
				G.p[u] = 0;
				term[u] = 0;
				en--;
				ok = 1;
				break;
			}
		}		
	}
	return ok;
}


int preprocessing(Graph &G){
	int n = G.V;
	int m = G.E;
	vi p = G.p;
	vvii adj = G.adj;

	vector<bool> term(n + 1);
	for (int v = 1; v <= n; v++)
		term[v] = G.p[v] > 0;

	vector< vector<int> > dist(n + 1, vector<int>(n + 1));

	int en = n; // excluded n, avoid conflicts insede the process since we do not remap
				// excluded verticies

	bool ok = 1;
	while (ok){
	
		do{
			do{
				do{
					do{
						ok = 0;
						ok |= pcdist(n, m, p, term, adj);
					} while (ok);
					
				ok |= least_cost(n, m, adj, dist);
				ok |= degree_1(n, m, adj, G);	

				} while(ok);
			
				// Now, we are going to add some stronger reductions that 	
				// need tranformations
				
				ok |= degree_2(n, m, en, term, adj);	

			} while (ok);	
		
			ok |= degree_3(n, m, en, term, adj, dist);
	
		} while (ok);

		// quebrado, não sei pq, olhar depois
		
	//	ok |= min_adj(n, m, en, term, adj, G);	
	}	
	
	int sum = 0;
	for (int v = 1; v <= G.V; v++)
		sum += adj[v].size();

	if (sum != check(G.V, adj))
		cout<<"Fudeu"<<endl;

	// cout<<"Vtx : "<<G.V<<' '<<en<<' '<<100*((double) (G.V - en)/G.V)<<"%"<<endl;
	// cout<<"Edg : "<<G.E<<' '<<m<<' '<<100*((double) (G.E - m)/G.E)<<"%"<<endl;



	G.adj = adj;
	G.E = m;
	//G.V = n; // é errado fazer isso sem modificar a indexação
	return en;	
=======
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

>>>>>>> c26894eb66ef944092aceb5c8be65419d1d78516
}

#ifndef PCST

int main(){
	string s = "../instances/PCSPG-JMP/P400.3.stp";
	//string s = "../instances/PCSPG-CRR/C02-A.stp";

<<<<<<< HEAD
	//string s = "../instances/morete.stp";

	Graph G;
	
	cout<<setprecision(2)<<fixed;

	if (STP_reader(s, G) != 0){
		cout<<"Arquivo não existe"<<endl;
		return 0;
	}

	preprocessing(G);
}

#endif
=======
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
>>>>>>> c26894eb66ef944092aceb5c8be65419d1d78516
