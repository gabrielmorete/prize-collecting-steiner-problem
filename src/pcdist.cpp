/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
   Computes prize-constrained distance
*/

#include"prepro.h"
#include"debug.h"
#include"stp_reader.h"
#include"graph.h"
using namespace std;

bool dijkstra(int src, int snk, int cst, vector<bool> &term, vi &p, vvii &adj, vi &dist){
	fill(all(dist), INF);	

	vector<bool> forbid(dist.size(), 0);

	priority_queue< pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>> > pq;

	pq.push({0, src});
	forbid[src] = true;

	int dst, vtx, w, c;

	while (!pq.empty()){
		tie(dst, vtx) = pq.top();
		pq.pop();
		
		if (dist[vtx] < dst)
			continue;
		
		if (term[vtx])
			forbid[vtx] = 1;

		dist[vtx] = dst;

		for (auto x : adj[vtx]){
			tie(w, c) = x;
			
			if (vtx == src and w == snk)
				continue;

			if (!forbid[w] and dst + c <= cst and dst + c - p[w] < dist[w]){
				if (w == snk){
					// dbg((vtx != src));
					// chapa;
					return true;
				}
				pq.push({ max(0, c + dst - p[w]), w});
			}
		}
	}

	return false;
}

bool pcdist(int &n, int &m, vi &p, vector<bool> &term, vvii &adj){
	bool ok = 0;

	vi dist(n + 1), dist2(n + 1);	

	int extra = 0;

	for (int v = 1; v <= n; v++)
		for (int i = 0; i < adj[v].size(); i++){
			if (dijkstra(v, adj[v][i].first, adj[v][i].second, term, p, adj, dist)){			
				erase_adj(v, adj[v][i].first, adj);

				swap(adj[v][i], adj[v][adj[v].size() - 1]);
				adj[v].pop_back();
				m--;
				ok = 1;
			}
			else if (dijkstra(adj[v][i].first, v, adj[v][i].second, term, p, adj, dist2)){
				erase_adj(v, adj[v][i].first, adj);

				swap(adj[v][i], adj[v][adj[v].size() - 1]);
				adj[v].pop_back();
				m--;
				ok = 1;
			}
			else{
				for (int u = 1; u <= n; u++)
					if (dist[u] < INF and dist2[u] < INF and dist[u] + dist2[u] + p[u] <= adj[v][i].second){
						erase_adj(v, adj[v][i].first, adj);

						swap(adj[v][i], adj[v][adj[v].size() - 1]);
						adj[v].pop_back();
						m--;
						extra++;
						ok = 1;
						break;
					}
			} 
		}			

	// dbg(extra);
	// dbg(check(n, adj));

	return ok;
}


// int main(){
// 	string s = "../instances/PCSPG-JMP/P400.3.stp";
// 	//string s = "../instances/PCSPG-CRR/C02-A.stp";

// 	//string s = "../instances/morete.stp";

// 	Graph G;
	
// 	cout<<setprecision(2)<<fixed;

// 	if (STP_reader(s, G) != 0){
// 		cout<<"Arquivo não existe"<<endl;
// 		return 0;
// 	}

// 	int n = G.V;
// 	int m = G.E;
// 	vvii adj = G.adj;
// 	vi p = G.p;

// 	vector<bool> term(n + 1, 0);
	
// 	for (int v = 1; v <= n; v++)
// 		term[v] = p[v] > 0;

// 	pcdist(n, m, p, term, adj);
// 	cout<<"Edg : "<<G.E<<' '<<m<<' '<<100*((double) (G.E - m)/G.E)<<"%"<<endl;

// 	while(G.E < m){
// 		m = G.E;
// 		pcdist(n, m, p, term, adj);
// 	}
// 	G.adj = adj;
// 	G.E = m;
	
// }