#ifndef MST
#define MST

#include <tuple>
#include "frac.h"
double reduction(frac &opt, vector<int> &v, vector<pair<int, int>> &edges, Graph &G){
	if (v.size() == 1)
		return 0;

	frac sum = frac(0);

	for (int u = 1; u <= G.V; u++)
		sum = sum + G.p[u];
	for (int u : v)
		sum = sum - G.p[u];
	
	vector<bool> in_sol(G.V + 1, 0);
	for (int u : v)
		in_sol[u] = 1;

	vector<tuple<frac, int, int>> edg;

	for (int u = 1; u <= G.V; u++)
		for (auto e : G.adj[u])
			if (in_sol[u] and in_sol[e.first])
				edg.push_back({e.second, u, e.first}); // only edges inside the vertex induced subgraph

	
	sort(edg.begin(), edg.end());
	
	edges.clear();

	Dsu dsu = Dsu(G.V + 1);

	frac new_opt = sum;
	for (auto e : edg){
		frac val;
		int w, x;
		tie(val, w, x) = e;

		if (dsu.merge(w, x)){
			new_opt = new_opt + val;
			edges.push_back({w, x});
		}
	}		

	double redu = ((opt - new_opt).val()/opt.val()) * 100;
	cout<<"reduced : "<<redu<<endl;

	opt = new_opt;
	return redu;
}

#endif