#ifndef DINIC
#define DINIC

#include <vector>
#include <cstdlib>	
#include <cstring>
#include <cmath>

const double EPS = 1e-4;
int sign(double x){ return (x > EPS) - (x < -EPS); }

typedef vector< vector<double> > matrix;	

const int nmax = 400;

const int MAXN = 6e3;
const int MAXM = 2e4;

int ned, first[MAXN], work[MAXN], dist[MAXN], q[MAXN]; // vertex information, integer

double cap[MAXM], bcap[MAXM]; // edge capacity
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
		if (dist[v] == dist[u] + 1 && sign(cap[e]) > 0){
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
		while (f = dfs(s, INF, t) and sign(f) > 0) 
			flow += f;
	}
	return flow;
}

#endif