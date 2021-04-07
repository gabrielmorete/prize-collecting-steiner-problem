#include<iostream>
#include<algorithm>
#include<vector>
#include<map>
#include<queue>

using namespace std;

struct frac {
	long long a, b;
	frac(int _a, int _b){ 
		int c = __gcd(__a, __b); 
		a = __a / c;
		b = __b / c;
	}

	frac operator+(frac p){
		long long top = a * p.b + b * p.a;
		long long bot = b * p.b;
		long long mdc = __gcd(top, bot);
		return {top / mdc, bot / mdc}; 
	}
	frac operator-(frac p){
		long long top = a * p.b - b * p.a;
		long long bot = b * p.b;
		if (top == 0) return {0, 1};
		long long mdc = __gcd(top, bot);
		return {top / mdc, bot / mdc}; 
	}
	frac operator*(frac p){
		long long top = a * p.a;
		long long bot = b * p.b;
		long long mdc = __gcd(top, bot);
		return {top / mdc, bot / mdc}; 
	}
	frac operator*(long long p){
		long long top = a * p;
		long long bot = b;
		long long mdc = __gcd(top, bot);
		return {top / mdc, bot / mdc}; 
	}
	frac operator/(long long p){
		long long top = a;
		long long bot = b * p;
		long long mdc = __gcd(top, bot);
		return {top / mdc, bot / mdc}; 
	}


};



const double INF = 1e100;

typedef priority_queue< pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>> > min_heap;

int n, N, mxatv;
vector<bool> u, l;
vector<int> o;
vector<double> dlt, d;
vector<vector<int>> L;
vector<pair<int, int>> edges;

map<pair<int, int>, double> cst;
map<pair<int, int>, pair<int, int>> A;

vector<set<double, int>> H0, H1;



inline double residual_cost(int v, int u){	return cst[{u, v}] - d[v] - d[u];}


inline double key(int i, int j){
	if (A[{i, j}].empty()) return INF;
	return residual_cost(A[{i,j}].first, A[{i,j}].second);
}


void init(){
	n = G.V;
	u.resize(n + 1);
	l.resize(n + 1);
	d.resize(n + 1);
	o.resize(n + 1);
	dlt.resize(n + 1);
	L.resize(n + 1);
	H.resize(n + 1);

	for (int v = 1; v <= n; v++){
		d[v] = 0;
		L[v].pb(v);
		o[v] = v;
		u[v] = l[v] = 1;;
		dlt[v] = p[v];
	}
	
	for (int i = 2; i <= n; i++)
		for (auto v : L[i])
			for (auto e : G.adj[v]){
				int u = e.first;
				int j = 0;
				if (o[u] == i)
					j = o[v];
				else
					j = o[u];
				if (key(i, j) > residual_cost(v, u))
					A[{i, j}] = A[{j, i}] = {u, v};		 
			}

	for (int i = 1; i <= n; i++)
		for (int j = 1; j < n; j++)
			if (j != i)
				if (!A[{i, j}].empty())
					H1[i].insert({key(i, j), j});		
}

void iterate(){
	double eps1, eps2;
	eps1 = eps2 = INF;

	int paux1 = 1, paux2 = 1, qaux2 = 1;
	for (int p = 1; p <= N; p++)
		if (u[p] == 1 and l[p] == 1){
			if (eps1 > dlt[p]){
				eps1 = dlt[p];
				paux1 = p;
			}

			if (!H0[p].empty()){
				if (eps2 > *H0[p].begin().first){
					tie(eps2, qaux2) = *H[0].begin()
					paux2 = p;
				}
			}

			if (!H1[p].empty()){
				if (2.0*eps2 > *H1[p].begin().first){
					eps2 = *H1[p].begin().first/2.0;
					qaux2 = *H1[p].begin().second;
					paux2 = p;
				}
			}
		}

	double eps = min(eps1, eps2);	

	for (int p = 1; p <= N; p++)
		if (u[p] and l[p] == 1){
			dlt[p] -= eps;
			for (int v : L[p])
				d[v] += eps
		}
	if (eps1 <= eps2)
		subcase1b(paux1);
	else
		subcase1a();		
}

void subcase1b(int p){
	l[p] = 0;
	mxatv--;
	for (int i = 1; i <= n; i++)
		if (i != p){
			if (H1[i].find({key(i, p), p}) != H1[i].end()){ // double problem, not good
				H1[i].erase({key(i, p), p});
				H0[i].insert({key(i, p), p});
			}
		}
}

void subcase1a(int p, int q){
	edges.push_back(A[{p, q}]);
	N++;
	L[N] = L[p];
	for (auto x : L[q])
		L[n].push_back(x);

	u[p] = u[q] = 0;
	u[N] = 1;

	if (l[q] == 1)
		mxatv--;

	dlt[N] = dlt[p] + dlt[q];

	l[N] = 1;
	for (int i = 1; i < N; i++)
		if (u[i] == 1){
			if (key(p, i) < key(q, i))
				A[{N, i}] = A[{i, N}] = A[{p, i}];
			else
				A[{N, i}] = A[{i, N}] = A[{q, i}];
		}
}


void prune(){

}

void pcst(){
	init();
	N = mxatv = n;
	while (mxatv > 1)
		iterate();
	prune();
}

int main(){

}