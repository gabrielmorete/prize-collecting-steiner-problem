/* Author : Gabriel Morete de Azevedo

   Event Handler for the PCST problem, Gutsfield
   algorithm for equivalent flow tree and Dinitz
   algorithm for maxflow.
*/

#include<vector>

const int MAXN = 1e5 + 3000;
const int MAXM = 9e5;
using namespace std;

const long double EPS = 1e-8;

namespace gmhu{

	int sign(coord x) { return (x > EPS) - (x < -EPS); }

	int ned, first[MAXN], work[MAXN], dist[MAXN], q[MAXN];
	int cap[MAXM], to[MAXM], nxt[MAXM];

	void init(){
	   memset(first, -1, sizeof first);
	   ned = 0;
	}

	void add(int u, int v, int f){
		to[ned] = v, cap[ned] = f;
		nxt[ned] = first[u];
		first[u] = ned++;
		
		to[ned] = u, cap[ned] = 0;
		nxt[ned] = first[v];
		first[v] = ned++;
	}

	int dfs(int u, int f, int t){
		if (u == t) return f;
		for (int &e = work[u]; e != -1; e = nxt[e]){
			int v = to[e];
			if (dist[v] == dist[u] + 1 && cap[e] > 0){
				int df = dfs(v, min(f, cap[e]), t);
				if (df > 0){
					cap[e] -= df;
					cap[e^1] += df;
					return df;
				}
			}
		}
		return 0;
	}

	bool bfs(int s, int t){
		memset(&dist, -1, sizeof dist);
		dist[s] = 0;
		int st = 0, en = 0;
		q[en++] = s;
		while (en > st){
			int u = q[st++];
			for (int e = first[u]; e != -1; e = nxt[e]){
				int v = to[e];
				if (dist[v] < 0 && cap[e] > 0){
					dist[v] = dist[u] + 1;
					q[en++] = v;
				}
			}
		}
		return dist[t] >= 0;
	}

	int dinic(int s, int t){
		int flow = 0, f;
		while (bfs(s, t)){
			memcpy(work, first, sizeof work);
			while (f = dfs(s, INF, t)) 
				flow += f;
		}
		return flow;
	}

	int n, m;
	int back_cap[MAXM]; // backup da capacidade
	int pai[MAXN], fpai[MAXN], mflow[MAXN][MAXN]; // tree
	int cuts[MAXN];

	void find_cut(int v){ // constroi o corte
		cuts[v] = true;
		for (int e = first[v]; e != -1; e = nxt[e])
			if (!cuts[to[e]] and cap[e] > 0)
				find_cut(to[e]);
	}

	int min_cut(int s, int t){
		memset(cuts, 0, sizeof cuts);
		memcpy(cap, back_cap, sizeof cap);

		int flow = dinic(s, t); // qualquer algoritmo de fluxo
		find_cut(s);			// que computa uma rede residual
								// basta acoplar o dinic desse caderno
		return flow;
	}

	void gomory_hu(){
		memset(mflow, -1, sizeof mflow); // -1 não tem corte que separa
		memcpy(back_cap, cap, sizeof cap);
		fill(pai, pai + n + 1, 1); // começa com k_{1, n - 1}

		for (int s = 2; s <= n; s++){
			int t = pai[s];
			fpai[s] = min_cut(s, t);
			
			mflow[s][t] = mflow[t][s] = fpai[s];

			for (int v = s + 1; v <= n; v++)
				if (cuts[v] and pai[v] == t)
					pai[v] = s;

			for (int v = 1; v < s; v++)
				if (mflow[s][v] == -1)
					mflow[s][v] = mflow[v][s] = min(mflow[s][t], mflow[t][v]);
		}
	}
}

double find_cut(int n, int r, double **adj, double *y){
	gmhu::init();
	vector<int> term, rterm(n + 1, 0);
	for (int v = 1; v <= n; v++)
		if (y[v] > 0.5){
			term.pb(v);
			rterm[v] = term.size() - 1;
		}
	for (int v = 1; v <= term.size(); v++)
		for (int u = 1; u <= term.size(); u++)
			if (u != v)
				gmhu::add(v, u, x[u][v]);
	gomory_hu();
}