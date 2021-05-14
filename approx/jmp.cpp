#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include "frac.h"
#include "graph.h"
#include "stp_reader.h"
#include "dsu.h"
using namespace std;

const int INF = 0x3f3f3f3f;
#define all(x) x.begin(), x.end()
#define dbg(x)  cout << #x << " = " << x << endl
#define chapa cout<<"oi meu chapa"<<endl

typedef pair<int, int> pii;

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

// Priority queue aux function

struct cmp{ // Comparador do set
	bool operator()(pii a, pii b) const {
		return key(a.first, a.second) < key(b.first, b.second);
	}
};

// O set vai ser de pair de int, o primeiro sendo o elemento, o segundo o conjunto
vector<set<pii, cmp>> H[2];

bool insert(int typ, int p, int q){
	H[typ][p].insert({q, p}); // elemento, conjunto
	return true;
}

bool remove(int typ, int p, int q){
	if (H[typ][p].count({q, p})){
		H[typ][p].erase({q, p});
		return true;
	}
	return false;
}

void copy(int typ, int p, int q){
	for (auto x : H[typ][p])
		H[typ][q].insert({x.first, q});
}

frac get_min_key(int typ, int p){
	int q = (*H[typ][p].begin()).first;
	return key(q, p);
}

int get_min_ele(int typ, int p){
	int q = (*H[typ][p].begin()).first;
	return q;
}

bool is_in(int typ, int p, int q){
	return H[typ][p].count({q, p}) > 0;
}


void subcase1b(int p);
void subcase1a(int p, int q);


void init(Graph G){
	n = G.V;
	u.resize(2 * n + 1);
	l.resize(2 * n + 1);
	d.resize(2 * n + 1);
	o.resize(2 * n + 1);
	dlt.resize(2 * n + 1);
	L.resize(2 * n + 1);
	H[0].resize(2 * n + 1);
	H[1].resize(2 * n + 1);

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
	
	for (int i = 2; i <= n; i++) // 2 mesmo?
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
		if (u[p] == 1 and l[p] == 1){ // atualizar isso aqui pode dar merda, corrigido
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
	// dbg(N);
	// dbg(mxatv);

	L[N] = L[p];
	for (auto x : L[q])
		L[N].insert(x);

	u[p] = u[q] = 0;
	u[N] = 1;

	if (l[q] == 1)
		mxatv--;

	dlt[N] = dlt[p] + dlt[q];

	l[N] = 1;
	for (int i = 1; i < N; i++) // errror here, fixed
		if (u[i] == 1){
			if (!A.count({p, i}) and !A.count({q, i})) // Not adjacent
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



vector< vector<int> > adj;
vector<frac> nw;

frac strong_prune(int v, int p){ // fix
	for (int u : adj[v]){
		if (u != p){
			frac f = strong_prune(u, v);
			if (cst[{v, u}] < f)
				nw[v] = nw[v] + nw[u] - cst[{v, u}];
		}
	}
	return nw[v];
}

frac opt;
vector<int> vertex;
void recover_prune(int v, int p){
	vertex.push_back(v);
	for (int u : adj[v]){
		if (u != p){
			if (cst[{v, u}] < nw[u]){
				edges.push_back({v, u});
				recover_prune(u, v);
			}
		}
	}
}

void check(Graph G){
	Dsu dsu = Dsu(n + 1);

	// for (auto e : edges)
	// 	cout<<e.first<<' '<<e.second<<endl;

	for (auto e : edges)
		if (!dsu.merge(e.first, e.second))
			cout<<"Circuito"<<endl;
	
	set<int> q;
	for (int i = 1; i <= n; i++)
		if (dsu.id[i] != i)
			q.insert(dsu.find(i));
	
	if (vertex.size() != edges.size() + 1)
		cout<<"Tamanho errado na árvore"<<endl;	

	if (q.size() > 1)
		cout<<"Mais de uma componente"<<endl;
	
	if (q.size() == 0 and vertex.size() > 1)
		cout<<"Vazio"<<endl;

	frac val = frac(0);
	for (int i = 0; i <= n; i++)
		val = val + G.p[i];

	for (int v : vertex)
		val = val - G.p[v];
	
	for (auto e : edges)
		val = val + cst[e];

	// val.print();
	// opt.print();

	if (val != opt)
		cout<<"Valor otimo distinto"<<endl;


	// val.print();
}


void prune(Graph G){
	nw.resize(n + 1);
	adj.resize(n + 1);

	for (auto e : edges){
		adj[e.first].push_back(e.second);
		adj[e.second].push_back(e.first);
	}

	frac sum = frac(0, 1);
	for (int i = 0; i <= n; i++){
		nw[i] = G.p[i];
		sum = sum + G.p[i];
	}

	for (int i = 1; i <= N; i++)
		if (l[i] == 1){ // only one
			int u = *L[i].begin();
			opt = sum - strong_prune(u, u);
			vertex.clear();
			edges.clear();
			recover_prune(u, u);
			break;
		}

	for (int v = 0; v <= n; v++){
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
	init(G);

	N = mxatv = n;
	while (mxatv > 1)
		iterate();
	
	// dbg(edges.size());

	prune(G);
	// dbg(vertex.size());
	// dbg(edges.size());
	check(G);
		

	return opt;
}

int main(){
	Graph G;
	// string fname = "../instances/PCSPG-JMP/P100.2.stp";
	string fname = "../instances/PCSPG-JMP/K100.stp";
	cout<<"Reader : "<<STP_reader(fname, G)<<endl;
	cout<<"OPT = "<<pcst(G)<<endl;
}