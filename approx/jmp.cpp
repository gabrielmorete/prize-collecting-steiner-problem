/* Author : Gabriel Morete de Azevedo
   Feofiloff at all implementation of the JMP algorithm
   Computes a 2-approximation of the PCST
*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include "frac.h"
#include "jmp.h"
#include "debug.h"
#include "graph.h"
#include "stp_reader.h"
#include "dsu.h"
#include "mst.h"
#include "preprocessing.cpp"
#include "tests.h"
using namespace std;

// const int INF = 0x3f3f3f3f;

int n, N, mxatv;
vector<bool> u, l;
vector<int> o;
vector<frac> dlt, d;

vector<set<int>> L; // Element sets
vector<pair<int, int>> edges;

map<pair<int, int>, frac> cst;
map<pair<int, int>, pair<int, int>> A;

inline frac residual_cost(int v, int u){ return frac(cst[{u, v}] - d[v] - d[u]); }

inline frac key(int i, int j){
	if (!A.count({i, j})) return frac(INF);
	return residual_cost(A[{i,j}].first, A[{i,j}].second);
}

struct cmp{ // Heap comparator, sort by heap id
	bool operator()(pair<int, int> a, pair<int, int> b) const {
		return key(a.first, a.second) < key(b.first, b.second);
	}
};

// Set constains a pair of ints, first is the element an second is the parent set L
vector<set<pair<int, int>, cmp>> H[2];
#include "heap.h" // heap auxiliary functions

void init(Graph G){
	n = G.V;
	u.resize(2 * n + 1, 0);
	l.resize(2 * n + 1, 0);
	d.resize(2 * n + 1, 0);
	o.resize(2 * n + 1, 0);
	dlt.resize(2 * n + 1, 0);
	L.resize(2 * n + 1);
	H[0].resize(2 * n + 1);
	H[1].resize(2 * n + 1);

	A.clear();
	cst.clear();
	edges.clear();
	for (int i = 0; i < 2 * n + 1; i++){
		L[i].clear();
		H[0][i].clear();
		H[1][i].clear();
	}

	for (int v = 1; v <= n; v++)
		for (auto e : G.adj[v])
			cst[{v, e.first}] = frac(e.second);

	for (int v = 1; v <= n; v++){
		d[v] = 0;
		L[v].insert(v);
		o[v] = v;
		u[v] = l[v] = 1;
		dlt[v] = G.p[v];
	}
	
	for (int i = 2; i <= n; i++)
		for (auto v : L[i])
			for (auto e : G.adj[v]){
				int u = e.first;
				int j = o[u];
				if (key(i, j) > residual_cost(v, u))
					A[{i, j}] = A[{j, i}] = {u, v};		 
			}

	for (int i = 1; i <= n; i++)
		for (int j = 1; j < n; j++)
			if (j != i and A.count({i, j}))
				insert(1, i, j);
}

void iterate(){
	frac eps1 = frac(INF);
	frac eps2 = frac(INF);
	
	int paux1 = 1, paux2 = 1, qaux2 = 1;
	for (int p = 1; p <= N; p++)
		if (u[p] == 1 and l[p] == 1){
			if (eps1 > dlt[p]){
				eps1 = dlt[p];
				paux1 = p;
			}

			if (!H[0][p].empty()){
				if (eps2 > get_min_key(0, p)){
					eps2 = get_min_key(0, p);
					qaux2 = get_min_ele(0, p);
					paux2 = p;
				}
			}	
	
			if (!H[1][p].empty()){
				if (eps2 * 2 > get_min_key(1, p)){
					eps2 = get_min_key(1, p);
					eps2 = eps2 / 2;
					qaux2 = get_min_ele(1, p);
					paux2 = p;
				}
			}
		}

	frac eps = min(eps1, eps2);	

	for (int p = 1; p <= N; p++)
		if (u[p] == 1 and l[p] == 1){
			dlt[p] = dlt[p] - eps;
			for (int v : L[p])
				d[v] = d[v] + eps;
		}

	if (eps1 <= eps2)
		subcase1b(paux1);
	else
		subcase1a(paux2, qaux2);		

}

void subcase1b(int p){
	l[p] = 0;
	mxatv--;
	for (int i = 1; i <= N; i++)
		if (i != p and remove(1, i, p))
			insert(0, i, p);
}

void subcase1a(int p, int q){
	edges.push_back(A[{p, q}]);
	N++;


	L[N] = L[p];
	for (auto x : L[q])
		L[N].insert(x);

	u[p] = u[q] = 0;
	u[N] = 1;

	if (l[q] == 1)
		mxatv--;

	dlt[N] = dlt[p] + dlt[q];

	l[N] = 1;
	for (int i = 1; i < N; i++)
		if (u[i] == 1){
			if (!A.count({p, i}) and !A.count({q, i}))
				continue;

			if (A.count({p, i}) and !A.count({q, i}))
				A[{N, i}] = A[{i, N}] = A[{p, i}];

			else if	(!A.count({p, i}) and A.count({q, i}))
				A[{N, i}] = A[{i, N}] = A[{q, i}];
			
			else if ((key(p, i) < key(q, i)))
				A[{N, i}] = A[{i, N}] = A[{p, i}];
			
			else
				A[{N, i}] = A[{i, N}] = A[{q, i}];
		}

	for (int h = 0; h <= 1; h++){ // New way
		vector<int> neigh;

		for (auto it : H[h][p])
			if (it.first != q)
				neigh.push_back(it.first);
		for (auto it : H[h][q])
			if (it.first != p)
				neigh.push_back(it.first);

		for (int u : neigh) // O(nlog(n))
			H[h][N].insert({u, N});
	}	

	for (int i = 1; i < N; i++) 
		if (u[i] == 1){
			remove(0, i, q);
			remove(1, i, q);
			remove(1, i, p);

			if (A.count({i, N}))
				insert(1, i, N);
		}
}

frac strong_prune(int v, int p, vector<frac> &nw, vector<vector<int>> &adj){
	for (int u : adj[v]){
		if (u != p){
			frac f = strong_prune(u, v, nw, adj);
			if (cst[{v, u}] < f)
				nw[v] = nw[v] + nw[u] - cst[{v, u}];
		}
	}
	return nw[v];
}

frac opt;
vector<int> vertex;
void recover_prune(int v, int p, vector<frac> &nw, vector<vector<int>> &adj){
	vertex.push_back(v);
	for (int u : adj[v]){
		if (u != p){
			if (cst[{v, u}] < nw[u]){
				edges.push_back({v, u});
				recover_prune(u, v, nw, adj);
			}
		}
	}
}

void check(Graph G){
	Dsu dsu = Dsu(n + 1);

	for (auto e : edges)
		if (!dsu.merge(e.first, e.second))
			cout << "\033[1;31mSolution is NOT a Tree : Cicle\033[0m\n";
	
	set<int> q;
	for (int i = 1; i <= n; i++)
		if (dsu.id[i] != i)
			q.insert(dsu.find(i));
	
	if (vertex.size() != edges.size() + 1)
		cout << "\033[1;31mSolution is NOT a Tree : |V| != |E| - 1\033[0m\n";

	if (q.size() > 1)
		cout << "\033[1;31mSolution is NOT a Tree : Disconnected\033[0m\n";
	
	if (q.size() == 0 and vertex.size() > 1)
		cout << "\033[1;31mChecking Error\033[0m\n";

	frac val = frac(0);
	for (int i = 0; i <= n; i++)
		val = val + G.p[i];

	for (int v : vertex)
		val = val - G.p[v];
	
	for (auto e : edges)
		val = val + cst[e];

	if (val != opt)
		cout << "\033[1;31mIncosistent optimal value\033[0m\n";
}


void prune(Graph G){
	vector< vector<int> > adj(n + 1);

	for (auto e : edges){
		adj[e.first].push_back(e.second);
		adj[e.second].push_back(e.first);
	}

	vector<frac> nw(n + 1);
	frac sum = frac(0, 1);
	for (int i = 1; i <= n; i++){
		nw[i] = G.p[i];
		sum = sum + G.p[i];
	}

	bool done = 0;
	for (int i = 1; i <= N; i++)
		if (l[i] == 1 and u[i] == 1){ // Should happen only once
			if (done)
				cout << "\033[1;31mMultiple active components!\033[0m\n";

			int u = *L[i].begin();
			opt = sum - strong_prune(u, u, nw, adj);
			vertex.clear();
			edges.clear();
			recover_prune(u, u, nw, adj);
			done = 1;
		}

	// Try a single vertex solution
	for (int v = 1; v <= n; v++){
		frac f = sum - G.p[v];
		if (opt > (sum - G.p[v])){
			opt = sum - G.p[v];
			vertex.clear();
			edges.clear();
			vertex.push_back(v);
		}
	}	
}


frac pcst(Graph G){
	// if (preprocessing(G) == -1)
	// 	cout << "\033[1;31mErro no preprocessamento\033[0m\n";

	init(G);



	N = mxatv = n;
	while (mxatv > 1)
		iterate();

	prune(G);

	check(G);

	// reduction(opt, vertex, edges, G);

	// check(G);

	return opt;
}

int main(){
	Graph G;
	// string fname = "../instances/PCSPG-JMP/P400.4.stp";
	// string fname = "../instances/PCSPG-JMP/K100.stp";
	// string fname = "../instances/ce.stp";

	// cout<<"Reader : "<<STP_reader(fname, G)<<endl;
	// STP_reader(fname, G);
	// cout<<"OPT = "<<pcst(G)<<endl;
	read_sol();
	runall();
	// // rune();

}