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


// Preprocessing function, if some error occour the return value is -1,
// othwewise is the number of real vertices remaining in the graph.
int preprocessing(Graph &G){
	int n = G.V;
	int m = G.E;
	vector<long long> p = G.p;
	vvii adj = G.adj;

	vector<bool> term(n + 1);
	for (int v = 1; v <= n; v++)
		term[v] = G.p[v] > 0;

	vector< vi > dist(n + 1, vi(n + 1));

	int en = n; // excluded n, avoid conflicts insede the process since we do not remap
				// excluded verticies

	bool ok = 1;
	// while (ok){
		ok = 0;
		// ok |= pcdist(n, m, p, term, adj);
		ok |= least_cost(n, m, adj, dist);
		ok |= degree_1(n, m, adj, G);	
	// // }	
	
	int sum = 0;
	for (int v = 1; v <= G.V; v++)
		sum += adj[v].size();

	if (sum != check(G.V, adj))
		return -1;
	
	// cout<<"Edg : "<<G.E<<' '<<m<<' '<<100*((double) (G.E - m)/G.E)<<"%"<<endl;

	G.adj = adj;
	G.p = p;
	G.E = m;
	return en;	
}

// #ifndef PCST

// int main(){
// 	// string s = "../instances/mst.stp";
// 	string s = "../instances/PCSPG-JMP/P400.4.stp";
// 	// string s = "../instances/E/E01-A.stp";

// 	Graph G;
	
// 	// cout<<setprecision(2)<<fixed;

// 	if (STP_reader(s, G) != 0){
// 		cout<<"Arquivo não existe"<<endl;
// 		return 0;
// 	}

// 	preprocessing(G);
// }

// #endif
