/* Author : Gabriel Morete de Azevedo

   Gutsfield algorithm for equivalent flow tree 
   and Dinitz algorithm for maxflow.
*/

#include<vector>
#include<cstring>
#include<iostream>
#include<cassert>
// #include "../src/stp_reader.cpp"
// #include "../src/debug.h"

// const int INF = 0x3f3f3f3f;

const int MAXN = 6e3;
const int MAXM = 2e4;
using namespace std;

const long double EPS = 1e-8;

int sign(double x) { return (x > EPS) - (x < -EPS); }

namespace gmhu{
	int ned, first[MAXN], work[MAXN], dist[MAXN], q[MAXN]; // vertex information, integer
	
	double cap[MAXM]; // edge capacity
	int to[MAXM], nxt[MAXM]; // edge information

	void init(){
	   memset(first, -1, sizeof first);
	   ned = 0;
	}

	void add(int u, int v, double f){
		to[ned] = v, cap[ned] = f;
		nxt[ned] = first[u];
		first[u] = ned++;
		
		to[ned] = u, cap[ned] = 0;
		nxt[ned] = first[v];
		first[v] = ned++;
	}

	double dfs(int u, double f, int t){
		if (u == t) 
			return f;
		for (int &e = work[u]; e != -1; e = nxt[e]){
			int v = to[e];
			if (dist[v] == dist[u] + 1 && cap[e] > 0){
				double df = dfs(v, min(f, cap[e]), t);
				if (sign(df) > 0){
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
				if (dist[v] < 0 and sign(cap[e]) > 0){
					dist[v] = dist[u] + 1;
					q[en++] = v;
				}
			}
		}

		return dist[t] >= 0;
	}

	double dinic(int s, int t){ // tested Auizu, ok!
		double flow = 0, f;
		while (bfs(s, t)){
			memcpy(work, first, sizeof work);
			while (f = dfs(s, INF, t)) 
				flow += f;
		}
		return flow;
	}

	int n, m;
	double back_cap[MAXM]; // backup da capacidade
	int pai[MAXN]; 
	double fpai[MAXN], mflow[MAXN][MAXN]; // tree
	int cuts[MAXN];

	void find_cut(int v){ // constroi o corte
		cuts[v] = true;
		for (int e = first[v]; e != -1; e = nxt[e])
			if (!cuts[to[e]] and sign(cap[e]) > 0)
				find_cut(to[e]);
	}

	double min_cut(int s, int t){
		memset(cuts, 0, sizeof cuts);
		memcpy(cap, back_cap, sizeof cap);

		int flow = dinic(s, t); // qualquer algoritmo de fluxo
		find_cut(s);			// que computa uma rede residual
								// basta acoplar o dinic desse caderno
		return flow;
	}

	void gomory_hu(){ // tested
		memcpy(back_cap, cap, sizeof cap);
		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
				mflow[i][j] = -1;

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


	int nterm, term[MAXN], rterm[MAXN]; // terminal and reverse permutation

	bool find_min_cut(int _n, double **adj, double *y){
		init();
		
		nterm = 0; // number of terminals, indexed by 1
		n = _n;	
		
		// frr(v, n)
		// 	dbg(y[v]);

		for (int v = 1; v <= n; v++)
			if (y[v] > 0.5){
				//dbg(v);
				term[++nterm] = v;
				rterm[v] = nterm;
			}


		for (int v = 1; v <= nterm; v++) // indexed by terminals
			for (int u = 1; u <= nterm; u++)
				if (u != v and sign(adj[term[u]][term[v]]) > 0){
	//				cout<<term[v]<<' '<<term[u]<<' '<<adj[term[v]][term[u]]<<endl;
					add(v, u, adj[term[v]][term[u]]);
				}
//		chapa;
		gomory_hu();

		// frr(v, nterm){
		// 	frr(u, nterm)
		// 		cout<<mflow[v][u]<<' ';
		// 	chapa;	
		// 	gnl;	
		// }


		// for (int v = 1; v <= nterm; v++){
		// 	cout<<term[v]<<' '<<v<<' '<<gmhu::pai[v]<<' '<<gmhu::fpai[v]<<endl;
		// }



		// now we want the global minimun cut

		int vtx1 = 0, vtx2 = 0;
		double minc = INF;
		for (int v = 1; v <= nterm; v++)
			for (int u = 1; u <= nterm; u++)
				if (u != v  and sign(mflow[v][u] - minc) < 0){
					minc = mflow[v][u];
					vtx1 = v;
					vtx2 = u; 
				}


		//cout<<"mincut : "<<vtx1<<' '<<vtx2<<' '<<mflow[vtx1][vtx2]<<endl;		

		if (vtx1 == 0){
			cout<<"Erro: no vertex was found"<<endl;
			return true;
		}

		if (sign(minc - 1) >= 0)
			return false; // found optimal solution

		min_cut(vtx1, vtx2);

		// cout<<"cuts : ";
		// for (int v = 1; v <= nterm; v++)
		// 	if (cuts[v] == 1)
		// 		cout<<term[v]<<' ';
		// gnl;	

		return true;
	}
}


// int main(){

// 	Graph G;

// 	string file_name = "ljubic.stp";
// 	if (int _code = STP_reader(file_name, G) != 0){
// 		cout<<"Error reading file - Code "<<_code<<endl;
// 		return 0;
// 	}

// 	int n = G.V, m = G.E;
// 	vector<vector<double> > adj(n + 1, vector<double>(n + 1));

// 	for (int v = 1; v <= n; v++)
// 		for (int u = 1; u <= n; u++)
// 			adj[v][u] = (u == v ? 0 : -1);

// 	for (int v = 1; v <= n; v++)
// 		for (auto u : G.adj[v])
// 			adj[v][u.first] = u.second;

// 	// frr(v, n){
// 	// 	frr(u, n)
// 	// 		cout<<adj[v][u]<<' ';
// 	// 	gnl;	
// 	// }

// 	double y[n + 1];

// 	fr(i, n + 1)
// 		y[i] = 1;
// 	y[1] = 0;	

// 	cout<<(gmhu::find_min_cut(n, 0, adj, y))<<endl;	
// }

// int main(){
// 	gmhu::init();
// 	cin>>gmhu::n>>gmhu::m;

// 	int x, y, z;

// 	fr(i, gmhu::m){
// 		cin>>x>>y;//>>z;
// 		gmhu::add(x, y, 1);
// 		gmhu::add(y, x, 1);
// 	}

// 	gmhu::gomory_hu();		

// 	int ans = 0;

// 	for (int i = 1; i <= gmhu::n; i++)
// 		for (int j = 1; j < i; j++)
// 			ans += gmhu::mflow[j][i];
// 	cout<<ans<<endl;	
// }

