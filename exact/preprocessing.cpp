/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

#include "prepro.h"
#include "stp_reader.h"
#include "graph.h"
#include "pcdist.cpp"
#include "dsu.h"
#include "debug.h"
#include <map>
#include <chrono>
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

bool mst_test(int &n, int &m, vector<bool> &term, vvii &adj, vi &p){
	bool tried[n + 1];	
	int pnew[n + 1], min_cst[n + 1];
	Dsu dsu = Dsu(n);

	vector< tuple<int, int, int> > edges;
	for (int v = 1; v <= n; v++)
		for (auto e : adj[v])
			if (v < e.first)
				edges.push_back({e.second, v, e.first}); 

	sort(edges.begin(), edges.end());		

	for (int v = 1; v <= n; v++)
		pnew[v] = p[v];

	bool okk = 1;
	while (okk){
		okk = 0;
		memset(tried, 0, sizeof tried);

		for (auto e : edges){
			int v, u, cs;
			tie(cs, v, u) = e;
			int rv = dsu.find(v);
			int ru = dsu.find(u);
			if (min(pnew[rv], pnew[ru]) > cs and ru != rv){
				if ((tried[rv] and cs != min_cst[rv]) or (tried[ru] and cs != min_cst[ru]))
					continue; // edge does not have minimal cost
				dsu.merge(u, v);
				pnew[dsu.find(v)] = pnew[rv] + pnew[ru] - cs;
				if (tried[rv] or tried[ru]){
					tried[dsu.find(rv)] = 1;
					min_cst[dsu.find(rv)] = min(tried[rv]? min_cst[rv] : INF, tried[ru]? min_cst[ru] : INF);
				}
				okk = 1;
			}
			else if (ru != rv){ // coudt add this edge, now it is minimal cost incident
				if (!tried[rv]){
					tried[rv] = 1;
					min_cst[rv] = cs;
				}
				if (!tried[ru]){
					tried[ru] = 1;
					min_cst[ru] = cs;
				}
			}
		}
	}	

	vector< tuple<int, int, int> > edges2;

	for (auto e : edges)
		if (dsu.find(get<1>(e)) != dsu.find(get<2>(e)))
			edges2.push_back(e);

	int removed = 0;	

	map< pair<int, int>, int > join;

	for (auto e : edges2){
		int v, u, cs;
		tie(cs, v, u) = e;	

		if (!join.count({dsu.find(v), dsu.find(u)})){
			join[{dsu.find(v), dsu.find(u)}] = 1;
			join[{dsu.find(u), dsu.find(v)}] = 1;
		}
		else {
			erase_adj(v, u, adj);
			erase_adj(u, v, adj);
			m--;
			removed++;
		}
	}	

	if(removed != 0)
		cout<<"mst : "<<removed<<endl;

	return removed > 0;
}

bool min_adj(int &n, int &m, int &en, vector<bool> &term, vvii &adj, vi &p){
	bool ok = 0;
	int cnt = 0, rm = 0;

	for (int v = 1; v <= n; v++){
		if (!term[v] or adj[v].empty())
			continue;

		int mn = get_min(v, adj);


		for (auto x : adj[v]){
			int u, cs;
			tie(u, cs) = x;
			
			if (min(p[v], p[u]) > cs and mn == cs){
				rm += fuse(v, u, n, adj);
				p[v] = p[v] + p[u] - cs;
				p[u] = 0;
				term[u] = 0;
				en--;
				cnt++;
				ok = 1;
				break;
			}
		}
	}
	// cout<<"sai do loop"<<endl;

	m -= rm;

	// cout<<"count   "<<cnt<<' '<<rm<<endl;

	return ok;
}


// Preprocessing function, if some error occour the return value is -1,
// othwewise is the number of real vertices remaining in the graph.
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
						// ok |= pcdist(n, m, p, term, adj);
					} while (ok);
					
				// ok |= least_cost(n, m, adj, dist);
				// ok |= degree_1(n, m, adj, G);	

				} while(ok);
			
				// Now, we are going to add some stronger reductions that 	
				// need tranformations
				
				// ok |= degree_2(n, m, en, term, adj);	

			} while (ok);	
		
			// ok |= degree_3(n, m, en, term, adj, dist);
	
		} while (ok);


		// ok |= min_adj(n, m, en, term, adj, p);


		ok |= mst_test(n, m, term, adj, p);	
	}	
	
	int sum = 0;
	for (int v = 1; v <= G.V; v++)
		sum += adj[v].size();

	if (sum != check(G.V, adj))
		return -1;
	
	// cout<<"Vtx : "<<G.V<<' '<<en<<' '<<100*((double) (G.V - en)/G.V)<<"%"<<endl;
	// cout<<"Edg : "<<G.E<<' '<<m<<' '<<100*((double) (G.E - m)/G.E)<<"%"<<endl;

	G.adj = adj;
	G.p = p;
	G.E = m;
	return en;	
}

#ifndef PCST

int main(){
	// string s = "../instances/mst.stp";
	// string s = "../instances/PCSPG-JMP/P400.4.stp";
	string s = "../instances/E/E01-A.stp";

	Graph G;
	
	cout<<setprecision(2)<<fixed;

	if (STP_reader(s, G) != 0){
		cout<<"Arquivo não existe"<<endl;
		return 0;
	}

	preprocessing(G);
}

#endif
