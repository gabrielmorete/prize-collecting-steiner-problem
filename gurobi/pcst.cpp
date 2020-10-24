/*	Author : Gabriel Morete de Azevedo
	Event Handler for the PCST

	Formulation:
	Given a network (G(V, E), c, p) with c = (c_e)_{e in E}, p = (p_v)_{v in V}
	define x \in E and y, r \in V binary such that for a tree solution T
		  / 1, if e \in E(T)	
	x_e = |
		  \ 0, otherwise

		  / 1, if v \in V(T)	
	y_v = |
		  \ 0, otherwise

		  / 1, if v is the root of T	
	r_v = |
		  \ 0, otherwise		  

	Hence, we have the following problem

		min cx + p(1 - y)
		st
		r*1 = 1 (only one root)
		y_v >= r_v (root is in the tree)
		x(\delta(W)) >= (1/|V|)y(W) -r(W), forwall W subset C
		x, y, r \in {0, 1}
*/

#include "gurobi_c++.h"
#include "../src/stp_reader.cpp"
#include "../src/debug.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }
void findsubtour(int n, double** sol, int* tourlenP, int* tour);

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.

class cut_tree: public GRBCallback{
	public:
		GBRVar* vary, varr;
		GRBVar** varx;
		int n;
		cut_tree(GRBVar *_y, GRBVar *_r, GRBVar** _x, int _n) {
			vary(_y), varr(_r), varx(_x), n(_n);
		}
	protected:
		void callback() {
			try {
				if (where == GRB_CB_MIPSOL) {
					// Found an integer feasible solution - does it visit every node?

					double y[n], r[n];
					double x[n + 1][n + 1];
					int i, j, len;
					for (int i = 1; i <= n; i++){
						y[i] = getSolution(y[i]);
						r[i] = getSolution(r[i]);
					}

					for (int i = 1; i <= n; i++)
						x[i] = getSolution(x[i], n);

					int root = 1, cnt = 0;
					for (int v = 1; v <= n; v++)
						if (r[v] > 0.5){
							root = r[v];
							cnt++;
						}
					assert(cnt != 1); // only one root	

					find_min_cut(n, root, x, y);

					if (len < n) {
						SUB_CNT++;//Estamos utilizando uma constraint de subtour elimination
						// Add subtour elimination constraint
						GRBLinExpr expr = 0;
						for (i = 0; i < len; i++)
							for (j = i+1; j < len; j++)
								expr += vars[tour[i]][tour[j]];
						addLazy(expr <= len-1);
					}

					for (i = 0; i < n; i++)
						delete[] x[i];
					delete[] x;
					delete[] tour;
					}
				} catch (GRBException e) {
					cout << "Error number: " << e.getErrorCode() << endl;
					cout << e.getMessage() << endl;
				} catch (...) {
					cout << "Error during callback" << endl;
				}
		}
};

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.  Result is returned in 'tour', and length is
// returned in 'tourlenP'.
	
void findsubtour(int n, double** sol, int* tourlenP, int* tour){
		bool* seen = new bool[n];
		int bestind, bestlen;
		int i, node, len, start;

		for (i = 0; i < n; i++)
			seen[i] = false;

		start = 0;
		bestlen = n+1;
		bestind = -1;
		node = 0;
		while (start < n) {
			for (node = 0; node < n; node++)
				if (!seen[node])
					break;
			if (node == n)
				break;
			for (len = 0; len < n; len++) {
				tour[start+len] = node;
				seen[node] = true;
				for (i = 0; i < n; i++) {
					if (sol[node][i] > 0.5 && !seen[i]) {
						node = i;
						break;
					}
				}
				if (i == n) {
					len++;
					if (len < bestlen) {
						bestlen = len;
						bestind = start;
					}
					start += len;
					break;
				}
			}
		}

	for (i = 0; i < bestlen; i++)
		tour[i] = tour[bestind+i];
	*tourlenP = bestlen;

	delete[] seen;
}

int main(){
	Graph G;

	string file_name = "ljubic.stp";


	if (int _code = STP_reader(file_name, G) != 0){
		cout<<"Error reading file - Code "<<_code<<endl;
		return 0;
	}

	int n = G.V, m = G.E;
	double adj[n + 1][n + 1];

	memset(adj, INF, sizeof adj);

	for (int v = 1; v <= n; v++)
		for (auto u : G.adj[v])
			adj[v][u.first] = u.second;

	for (int i = 1; i <= G.V; i++)
		for (pair<int, int> u : G.adj[i])
			if (u.first > i)
				cout<<u.first<<' '<<i<<endl;

	GRBEnv env = GRBEnv();
	GRBVar x[n + 1][n + 1];
	GRBVar y[n + 1];
	GRBVar r[n  + 1];

	
	try {

		GRBModel model = GRBModel(env);

		// We will use an event Handler

		model.set(GRB_IntParam_LazyConstraints, 1);

		// Create binary decision variables

		for (int i = 1; i < n; i++) {
			for (int j = 1; j <= i; j++) {
				vars[i][j] = model.addVar(0.0, (adj[i][j] != INF)? 1.0 : 0.0, adj[i][j], GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
				vars[j][i] = vars[i][j];
			}
		}

		GRBLinExpr expr = 0;
		for (int v = 1; v <= n; v++){
			y[v] = model.addVar(0.0, 1.0, -G.p[v], GRB_BINARY, "y_"+itos(v));
			y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_"+itos(v));
		}

		int sum = 0;
		for (int v = 1; v <= n; v++)
			sum += G.p[v];
		y[0] = model.addVar(1.0, 1.0, sum, GRB_BINARY, "sum_of_values")

		// Forbid edge from node back to itself
		for (int i = 1; i < n; i++)
			vars[i][i].set(GRB_DoubleAttr_UB, 0);

		// Root is in the tree
		for (int v = 1; v <= n; v++)
			model.addConstr(y[v] >= r[v], "in_tree_root_"+itos(v));

		// Only one root
		expr = 0;
		for (int v = 1; v <= n; v++)	
			expr += r[v];
		model.addConstr(expr == 1, "one_root");

		// Set callback function

		cut_tree cb = cut_tree(y, r, x, n);
		model.setCallback(&cb);

		// Optimize model

		model.optimize();

		// Extract solution

		if (model.get(GRB_IntAttr_SolCount) > 0) { //////////////////////////
			double **sol = new double*[n];
			for (int i = 0; i < n; i++)
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

			int* tour = new int[n];
			int len;

			findsubtour(n, sol, &len, tour);
			assert(len == n);

			cout << "Tour: ";
			for (int i = 0; i < len; i++)
				cout << tour[i] << " ";
			cout << endl;
			cout << "SUB_ELIM:"<<SUB_CNT<<endl;
			for (int i = 0; i < n; i++)
				delete[] sol[i];
			delete[] sol;
			delete[] tour;
		}

	} catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Error during optimization" << endl;
	}

	for (int i = 0; i < n; i++)
	delete[] vars[i];
	delete[] vars;
	delete env;
	return 0;
}
