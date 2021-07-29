/*	Author : Gabriel Morete de Azevedo
	Event Handler for the PCST

	Formulation:
	Given a network (G(V, E), c, p) with c = (c_e)_{e in E}, p = (p_v)_{v in V}
	define x \in E and y \in V binary such that for a tree solution T
		  / 1, if e \in E(T)	
	x_e = |
		  \ 0, otherwise

		  / 1, if v \in V(T)	
	y_v = |
		  \ 0, otherwise

	We add a new vertex r (v_0) to act as root of the tree, hence, we will solve
	the rooted version of the problem.

	Hence, we have the following formulation

		min cx + p(1 - y)
		
		subject to		
		
		y_r = 1 (root is in the tree)
		sum{v \in V} x_{rv} = 1 (only one neighbour for root)
		\sum_{e \in delta(v)} x_e <= |v|y_v, forall v \in V	
		x(\delta(W)) >= (1/|W|)y(W), forwall W subset V, r \in W (connectivity)
		x_{rv} <= 1 - y_u, forall u < v and v \in V (symmetry breaking)
		x, y \in {0, 1}

	Here we will try to use a simpler event handler (a dfs from the root) to
	find infeaseble cuts. Hence we will define the edges as integer and let
	Gurobi handle the integrality;

*/

#include "gurobi_c++.h"
#include "stp_reader.cpp"
#include "debug.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

const double EPS = 1e-6;
int sign(double x){ return (x > EPS) - (x < -EPS); }

typedef vector< vector<double> > matrix;	

// // returns the # of visited nodes
// int dfs(int v, int n, int *vis, double *y, matrix adj){
// 	vis[v] = 1;

// 	int tree = 1;
// 	for (int u = 1; u <= n; u++)
// 		if (y[u] > 0.5 and adj[v][u] > 0.5 and !vis[u])
// 			tree += dfs(u, n, vis, y, adj);

// 	return tree;	
// /}

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
	to[ned] = v, bcap[ned] = f;
	nxt[ned] = first[u];
	first[u] = ned++;
	
	to[ned] = u, bcap[ned] = 0;
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
				if (sign(cap[e] - df) <= 0){
					cap[e] = 0;
					cap[e^1] = bcap[e];
				}
				else{
					cap[e] -= df;
					cap[e^1] += df;
				}
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
	memcpy(cap, bcap, sizeof cap);
	double flow = 0, f;
	//cout<<"entrei"<<endl;
	while (bfs(s, t)){
		memcpy(work, first, sizeof work);
		while (f = dfs(s, INF, t) and sign(f) > 0) 
			flow += f;
	}
	//cout<<"sai"<<endl;
	return flow;
}


class cut_tree: public GRBCallback{
	public:
		GRBVar* vary; 
		GRBVar** varx;
		vector<vector<bool>> adj;
		int n;
		cut_tree(GRBVar* _y, GRBVar** _x, int _n, vector<vector<bool>> _adj) {
			vary =_y; varx = _x;  n = _n; adj = _adj;
		}
	protected:
		void callback() {
			try {
				if (where == GRB_CB_MIPSOL) {
					// Found an integer feasible solution - does it form a tree?

					double y[n + 1];
					vector<vector<double>> x(n + 1, vector<double>(n + 1));
					
					for (int v = 0; v <= n; v++)
						y[v] = getSolution(vary[v]);

					for (int v = 0; v <= n; v++){
						for (int u = 0; u <= n; u++)
							x[v][u] = getSolution(varx[v][u]);
					}

					init(); // start dinic	

					for (int v = 0; v <= n; v++)
						for (int u = 0; u <= n; u++)
							if (x[v][u] > 0){
								add(v, u, x[v][u]);
							}

					for (int vtx = 1; vtx <= n; vtx++){
						if (sign(dinic(0, vtx) - y[vtx]) >= 0)
							continue;

						GRBLinExpr cut = 0;

						for (int v = 0; v <= n; v++){
							if (dist[v] >= 0){ // terminal inside the cut
								for (int u = 1; u <= n; u++){
									if (dist[u] < 0)
										cut += varx[v][u];
								}
							}	
						}

						addLazy(cut >= vary[vtx]);



						// for (int v = 1; v <= n; v++)
						// 	if (vis[v] == 0 and y[v] > 0.5)
						// 		addLazy(cut >= vary[v]);

						// for (int v = 0; v <= n; v++){
						// 	if (vis[v] == 1){ // terminal inside the cut
						// 		for (int u = 1; u <= n; u++){
						// 			if (adj[v][u])
						// 				x[v][u] = 1;
						// 		}
						// 	}	
						// }
						// memset(vis, 0, sizeof vis);						
					}
				}
				} catch (GRBException e) {
					cout << "Error number: " << e.getErrorCode() << endl;
					cout << e.getMessage() << endl;
				} catch (...) {
					cout << "Error during callback" << endl;
				}
		}
};
using dtype = long long int;

int main(){
	Graph G; // type graph from stp reader

	//string file_name = "../instances/lujibic/ljubic3.stp";
	//string file_name = "../instances/B/b01.stp";
	//string file_name = "../instances/PCSPG-JMP/K100.10.stp";
	string file_name = "../instances/PCSPG-JMP/K400.stp";
	//string file_name = "../instances/PCSPG-JMP/P100.4.stp";

	int _code = STP_reader(file_name, G);
	if (_code != 0){
		cout<<"Error reading file - Code "<<_code<<endl;
		return 0;
	}

	// Add agressive preprocessing step
	////////////////////////////////////
	////////////////////////////////////


	int n = G.V, m = G.E;
	vector<vector<bool>> adj(n + 1, vector<bool>(n + 1, 0));
	bool term[n + 1];
	dtype cost[n + 1][n + 1], dist[n + 1][n + 1];

	G.p[0] = 0;

	term[0] = 1;
	for (int v = 1; v <= n; v++){
		term[v] = 0;
		if (G.p[v] > 0)
			term[v] = 1;
	}

	for (int v = 0; v <= n; v++)
		for (int u = 0; u <= n; u++){
			adj[v][u] = false;
			cost[v][u] = 0;
		}

	for (int v = 0; v <= n; v++){	
		adj[v][v] = true;
		adj[0][v] = true;
		cost[v][v] = 0;
		if (v > 0)
			cost[0][v] = -G.p[v];
	}	

	for (int v = 1; v <= n; v++) // graph is 1-indexed
		for (auto u : G.adj[v]){  // root is node 0
			adj[v][u.first] = true;
			cost[v][u.first] = u.second - G.p[u.first];
		}	

	GRBEnv env = GRBEnv();
	GRBVar **x = new GRBVar*[n + 1];
	for (int i = 0; i <= n; i++)
		x[i] = new GRBVar[n + 1];

	GRBVar *y = new GRBVar[n + 1];
	GRBVar *r = new GRBVar[n + 1];

	try {

		GRBModel model = GRBModel(env);

		// We will use an event Handler
		// Setting some parameters
		model.set(GRB_IntParam_LazyConstraints, 1);
		model.set(GRB_DoubleParam_TimeLimit, 1000.0);
		model.set(GRB_DoubleParam_Heuristics, 0.15);
		

		// Create binary decision variables
		for (int v = 0; v <= n; v++)
			for (int u = 0; u <= n; u++)
				x[v][u] = model.addVar(0.0, (double) adj[v][u], 0, GRB_BINARY, "x_"+itos(v)+"_"+itos(u));
			
		for (int v = 0; v <= n; v++)
			y[v] = model.addVar(0.0, 1.0, 0, GRB_BINARY, "y_"+itos(v));
		
		// Forbid edge from node back to itself
		for (int v = 0; v <= n; v++)
			x[v][v].set(GRB_DoubleAttr_UB, 0);

		// Forbid edge from node to the root
		for (int v = 0; v <= n; v++)
			x[v][0].set(GRB_DoubleAttr_UB, 0);

		// Add non-costumer verticies restiction
		for (int v = 1; v <= n; v++)
			if (!term[v]) // non costumer can't be root
				x[0][v].set(GRB_DoubleAttr_UB, 0);		

		// Root is in the tree
		y[0].set(GRB_DoubleAttr_LB, 1);

		// Only one root
		GRBLinExpr expr = 0;
		for (int v = 1; v <= n; v++)	
			expr += x[0][v];
		model.addConstr(expr == 1, "one_root");

		// Directed tree equivalence
		for (int v = 1; v <= n; v++){
			expr = 0;
			for (int u = 0; u <= n; u++)
				expr += x[u][v];
			model.addConstr(expr == y[v], "direct+"+itos(v));
		}

		// // Incident inside the tree
		// for (int v = 1; v <= n; v++){
		// 	expr = 0;
		// 	for (int u = 0; u <= n; u++)
		// 		expr += x[v][u];
		// 	model.addConstr((n + 1) * y[v] >= expr, "in_edge_in_tree"+itos(v));
		// }

		// Symmetry break
		for (int v = 1; v <= n; v++) // fix after non-costumer-restriction
			for (int u = 1; u < v; u++)
				if (term[u]) // non-costumer might mess
					model.addConstr(x[0][v] <= 1 - y[u], "symmetry_"+itos(v)+"_"+itos(u));

		// Only one orientation per edge		
		for (int v = 1; v <= n; v++)
			for (int u = 1; u <= n; u++)
				if (adj[v][u])
					model.addConstr(x[v][u] + x[u][v] <= y[v], "orientation_"+itos(v)+"_"+itos(u));

		// Flow-balance constraints
		for (int v = 1; v <= n; v++)
			if (!term[v]){
				expr = 0;
				for (int u = 1; u <= n; u++)
					expr += x[v][u] - x[u][v];
				model.addConstr(expr >= 0, "flow_balance"+itos(v));
			}		


		// Tree size
		// expr = 0;
		// for (int v = 1; v <= n; v++)
		// 	expr += y[v];
		// for (int v = 1; v <= n; v++)
		// 	for (int u = 1; u <= n; u++)
		// 		expr -= x[v][u];
		// model.addConstr(expr == 1, "tree_size");					

		// // Non-leaf non-costumer (Lucena)
		// for (int v = 1; v <= n; v++){
		// 	expr = 0;
		// 	for (int u = 1; u <= n; u++)
		// 		expr += x[v][u];
		// 	if (!term[v])				
		// 		model.addConstr(expr >= 2 * y[v]);
		// 	// else // Non-single vertex
		// 	// 	model.addConstr(expr >= y[v]);
		// }	

		// Preprocessing
		int pren, prem;
		pren = prem = 0;

		// for (int v = 1; v <= n; v++) // dist for floyd-warshall
		// 	for (int u = 1; u <= n; u++)
		// 		dist[v][u] = (adj[v][u] == inf? INF : adj[v][u]);

		// for (int k = 1; k <= n; k++) // floyd
		// 	for (int v = 1; v <= n; v++)
		// 		for (int u = 1; u <= n; u++)
		// 			dist[v][u] = min(dist[v][u], dist[v][k] + dist[k][u]);

		// for (int v = 1; v <= n; v++) // least cost test
		// 	for (int u = 1; u < v; u++)
		// 		if (adj[v][u] > dist[v][u]){
		// 			prem++;
		// 			adj[v][u] = adj[u][v] = inf;
		// 		}


		// for (int v = 1; v <= n; v++) // degree 1 test
		// 	if (G.adj[v].size() == 1 and (G.adj[v][0].second - G.p[v]) >= 0){
		// 		y[v].set(GRB_DoubleAttr_UB, 0);
		// 		x[v][G.adj[v][0].first].set(GRB_DoubleAttr_UB, 0);
		// 		pren++;
		// 		prem++;
		// 	}	

		// int a, b;
		// dtype cst;
		// for (int v = 1; v <= n; v++) // degree 2 test
		// 	if (!term[v] and G.adj[v].size() == 2){
		// 		a = G.adj[v][0].first;
		// 		b = G.adj[v][1].first;	
		// 		cst = adj[v][a] + adj[v][b];

		// 		pren++;
		// 		prem++;
		// 		if (adj[v][a] == inf or adj[v][b] == inf){
		// 			x[v][a].set(GRB_DoubleAttr_UB, 0);
		// 			x[v][b].set(GRB_DoubleAttr_UB, 0);
		// 			y[v].set(GRB_DoubleAttr_UB, 0);
		// 			prem++;
		// 		}
		// 		else if (adj[a][b] == inf or adj[a][b] > cst){
		// 			adj[a][b] = cst;
		// 			model.addConstr(x[a][b] == y[v]);
		// 		}
		// 		else{
		// 			prem++;
		// 			y[v].set(GRB_DoubleAttr_UB, 0);
		// 			x[v][a].set(GRB_DoubleAttr_UB, 0);
		// 			x[v][b].set(GRB_DoubleAttr_UB, 0);
		// 		}

		// 		adj[v][a] = adj[a][v] = inf;
		// 		adj[v][b] = adj[b][v] = inf;
		// 	}		

		cout<<"preprocessing : "<<pren<<" "<<prem<<endl;		

		// Exclude inexisting edges
		// for (int v = 1; v <= n; v++)
		// 	for (int u = 1; u <= v; u++)
		// 		if (!adj[v][u])
		// 			x[v][u].set(GRB_DoubleAttr_UB, 0);



		// Set Obj
		dtype obj = 0;
		for (int v = 1; v <= n; v++)
			obj += G.p[v];
		dbg(obj);
		expr = 0;
		for (int v = 0; v <= n; v++)
			for (int u = 0; u <= n; u++)
				if (adj[v][u]) 
					expr += cost[v][u] * x[v][u];

		expr +=	obj;

		model.setObjective(expr, GRB_MINIMIZE);	

		// Set callback function
		cut_tree cb = cut_tree(y, x, n, adj);
		model.setCallback(&cb);

		// Optimize model

		model.optimize();

		// Extract solution

		if (model.get(GRB_IntAttr_SolCount) > 0) { /////////////////////////
			double *sol = model.get(GRB_DoubleAttr_X, y, n + 1);
			
			double **edge = new double*[n + 1];
			for (int v = 0; v <= n; v++)
				edge[v] = model.get(GRB_DoubleAttr_X, x[v], n + 1);

			cout << "Tree: ";
			for (int i = 0; i <= n; i++)
				if(sol[i] >= 0.5)
					cout<<i<<' ';
			cout << endl;
			cout<<"Edges"<<endl;
			for (int v = 0; v <= n; v++){
				for (int u = 0; u <= n; u++){
					 if (edge[u][v] > 0.5)
						cout<<u<<' '<<v<<' '<<cost[u][v]<<endl;
				}
			}	
			cout<<"preprocessing : "<<pren<<" "<<prem<<endl;		


			for (int v = 0; v <= n; v++)
				delete[] edge[v];
			delete[] edge;			
			delete[] sol;
		}

	} catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Error during optimization" << endl;
	}

	for (int i = 0; i <= n; i++)
		delete[] x[i];
	delete[] x;
	delete[] y;
	delete[] r;

	return 0;
}
