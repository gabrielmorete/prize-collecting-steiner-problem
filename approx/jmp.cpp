/* Author : Gabriel Morete de Azevedo
   Feofiloff et al implementation of the JMP algorithm
   Computes a 2-approximation of the Prize-collecting
   steiner tree problem
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
#include "graph.h"
#include "stp_reader.h"
#include "dsu.h"
#include "mst.h"
#include "tests.h"
using namespace std;

#define dbg(x)  cout << #x << " = " << x << endl
#define endl '\n'
#define INF 0x3f3f3f3f

int n, N, mxatv;
vector<bool> u, l; // u[i] = i is maximal? l[i] = i is active?
vector<frac> dlt, d; // Delta value for set, d value for element

// Solution
frac opt;
vector<int> vertex;
vector<pair<int, int>> edges;

vector<vector<int>> L; // Element sets

map<pair<int, int>, frac> cst;
map<pair<int, int>, pair<int, int>> A;

inline frac residual_cost(int v, int u){ return cst[{u, v}] - d[v] - d[u]; }

inline frac key(int i, int j){
	if (!A.count({i, j})) return frac(INF);
	return residual_cost(A[{i,j}].first, A[{i,j}].second);
}

struct cmp{ // Heap comparator, sort by heap id
	bool operator()(pair<int, int> a, pair<int, int> b) const {
		if (key(a.first, a.second) == key(b.first, b.second)){
			return a.first < b.first;
		}
		return key(a.first, a.second) < key(b.first, b.second);
	}
};

// Set constains a pair of ints, first is the element an second is the parent set L
vector<set<pair<int, int>, cmp>> H[2];
#include "heap.h" // heap auxiliary functions

void init(Graph G){
	N = mxatv = n = G.V;
	u.resize(2 * n + 1, 0);
	l.resize(2 * n + 1, 0);
	d.resize(2 * n + 1, 0);
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
		L[v].push_back(v);
		u[v] = l[v] = 1;
		dlt[v] = frac(G.p[v]);
	}
	
	// Optimizing for instances with few terminals	
	for (int v = 1; v <= n; v++)
		if (G.p[v] == 0){
			l[v] = 0;	
			mxatv--;
		}

	for (int v = 1; v <= n; v++)
		for (auto e : G.adj[v]){
			int u = e.first;
			A[{v, u}] = A[{u, v}] = {u, v};		 
			if (l[u])
				insert(1, v, e.first);	
			else	
				insert(0, v, e.first);	
		}
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
				if ((eps2 * 2) > get_min_key(1, p)){
					eps2 = get_min_key(1, p);
					eps2 = eps2 / 2;
					qaux2 = get_min_ele(1, p);
					paux2 = p;
				}
			}
		}

	frac eps = min(eps1, eps2);	

	if (eps >= frac(INF)){
		cout << "\033[1;31mRuntime error. Invalid epsilon. Aborting!\033[0m\n";
		exit(-1);
	}

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

	L[N].resize(L[p].size() + L[q].size());
	int a = 0, b = 0, cnt = 0;
	while (a < L[p].size() and b < L[q].size()){
		if (L[p][a] <= L[q][b])
			L[N][cnt++] = L[p][a++];
		else
			L[N][cnt++] = L[q][b++];
	}

	while (a < L[p].size() and b == L[q].size())
			L[N][cnt++] = L[p][a++];
	
	while (a == L[p].size() and b < L[q].size())
			L[N][cnt++] = L[q][b++];

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

	for (int h = 0; h <= 1; h++){ // New way, O(nlog(n))
		for (auto it : H[h][p])
			if (it.first != q)
				insert(h, N, it.first);
			
		for (auto it : H[h][q])
			if (it.first != p)
				insert(h, N, it.first);
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

void check(Graph G){
	Dsu dsu = Dsu(n + 1);
	
	if (vertex.size() != edges.size() + 1)
		cout << "\033[1;31mSolution is NOT a Tree : |V| != |E| - 1\033[0m\n";

	for (auto e : edges)
		if (!dsu.merge(e.first, e.second))
			cout << "\033[1;31mSolution is NOT a Tree : Cicle\033[0m\n";
	
	set<int> q;
	for (int i = 1; i <= n; i++)
		if (dsu.id[i] != i)
			q.insert(dsu.find(i));

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
		cout << "\033[1;31mInconsistent optimal value\033[0m\n";
}

frac strong_prune(int v, int p, vector<frac> &nw, vector<long long> &pi, vector<vector<int>> &adj){
	nw[v] = frac(pi[v]);
	for (int u : adj[v]){
		if (u != p){
			frac f = strong_prune(u, v, nw, pi, adj);
			if (cst[{v, u}] < f)
				nw[v] = nw[v] + nw[u] - cst[{v, u}];
		}
	}
	return nw[v];
}

void recover_prune(int v, int p, vector<frac> &nw, vector<vector<int>> &adj){
	vertex.push_back(v);
	for (int u : adj[v]){
		if (u != p and cst[{v, u}] < nw[u]){
			edges.push_back({v, u});
			recover_prune(u, v, nw, adj);
		}
	}
}

void prune(Graph G){
	vector< vector<int> > adj(n + 1);

	for (auto e : edges){
		adj[e.first].push_back(e.second);
		adj[e.second].push_back(e.first);
	}

	vector<frac> nw(n + 1);
	frac sum = frac(0);
	for (int i = 1; i <= n; i++)
		sum = sum + G.p[i];

	opt = frac(INF);
	int id = 1;

	for (int v = 1; v <= n; v++){
		frac aux = sum - strong_prune(v, v, nw, G.p, adj);
		if (aux < opt){
			opt = aux;
			id = v;
		}
	}

	opt = sum - strong_prune(id, id, nw, G.p, adj);
	vertex.clear();
	edges.clear();
	recover_prune(id, id, nw, adj);
}


frac pcst(Graph G){
	init(G);

	while (mxatv > 1)
		iterate();

	prune(G);

	check(G);

	mst_reduction(opt, vertex, edges, G);

	check(G);

	return opt;
}

int main(){
	read_sol();
	runall();
}
