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
#include "min_cut.cpp"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.

class cut_tree: public GRBCallback{
	public:
		GRBVar* vary; 
		GRBVar* varr;
		GRBVar** varx;
		int n;
		cut_tree(GRBVar* _y, GRBVar* _r, GRBVar** _x, int _n) {
			vary =_y; varr = _r; varx = _x;  n = _n;
		}
	protected:
		void callback() {
			try {
				if (where == GRB_CB_MIPSOL) {
					// Found an integer feasible solution - does it visit every node?

					double y[n + 1], r[n + 1];
					double *x[n + 1];
					
					for (int i = 1; i <= n; i++){
						y[i] = getSolution(vary[i]);
						r[i] = getSolution(varr[i]);
					}

					for (int i = 1; i <= n; i++)
						x[i] = getSolution(varx[i], n + 1);

					int root = 1, cnt = 0;
					for (int v = 2; v <= n; v++)
						if (r[v] > r[root])
							root = r[v];

					for (int v = 1; v <= n; v++)
						cout<<y[v]<<' '<<r[v]<<endl;	

					for (int v = 0; v <= n; v++){
						for (int u = 0; u <= n; u++)
							cout<<x[u][v]<<' ';
						gnl;
					}



				
					if (gmhu::find_min_cut(n, x, y)){
						GRBLinExpr lhs = 0, rhs = 0;
						
						bool ok = !(gmhu::cuts[gmhu::rterm[root]]); // other side of the cut

						for (int v = 1; v <= gmhu::nterm; v++)
							for (int u = v + 1; u <= gmhu::nterm; u++)
								if (gmhu::cuts[v] != gmhu::cuts[u])
									lhs += varx[gmhu::term[v]][gmhu::term[u]]; // the cut

						double divi = 0;		
						for (int v = 1; v <= gmhu::nterm; v++)
							if (gmhu::cuts[v] == ok){ // side with no root
								divi += 1;	
							}	

						for (int v = 1; v <= gmhu::nterm; v++)
							if (gmhu::cuts[v] == ok){ // side with no root
								rhs += (vary[gmhu::term[v]]) / divi;
						}	
						for (int v = 1; v <= n; v++)
							rhs -= r[v];
						addLazy(lhs >= rhs);
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

int main(){
	Graph G; // type graph from stp reader

	string file_name = "ljubic.stp";

	int _code = STP_reader(file_name, G);

	if (_code != 0){
		cout<<"Error reading file - Code "<<_code<<endl;
		return 0;
	}

	int n = G.V, m = G.E;
	double adj[n + 1][n + 1];

	for (int v = 0; v <= n; v++)
		for (int u = 0; u <= n; u++)
			adj[v][u] = INF;

	for (int v = 1; v <= n; v++)
		for (auto u : G.adj[v])
			adj[v][u.first] = u.second;

	// for (int i = 1; i <= G.V; i++)
	// 	for (pair<int, int> u : G.adj[i])
	// 		if (u.first > i)
	// 			cout<<u.first<<' '<<i<<endl;

	GRBEnv env = GRBEnv();
	GRBVar **x = new GRBVar*[n + 1];
	for (int i = 0; i <= n; i++)
		x[i] = new GRBVar[n + 1];

	// 	vars = new GRBVar*[n];
	// for (i = 0; i < n; i++)
	// 	vars[i] = new GRBVar[n];

	GRBVar *y = new GRBVar[n + 1];
	GRBVar *r = new GRBVar[n + 1];

	try {

		GRBModel model = GRBModel(env);

		// We will use an event Handler

		model.set(GRB_IntParam_LazyConstraints, 1);

		// Create binary decision variables

		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= i; j++) {
				x[i][j] = model.addVar(0.0, 1.0, adj[i][j], GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
				x[j][i] = x[i][j];
			}
		}

		GRBLinExpr expr = 0;
		for (int v = 0; v <= n; v++){
			y[v] = model.addVar(0.0, 1.0, -G.p[v], GRB_BINARY, "y_"+itos(v));
			r[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_"+itos(v));
		}

		int sum = 0;
		for (int v = 1; v <= n; v++)
			sum += G.p[v];
		y[0] = model.addVar(1.0, 1.0, sum, GRB_BINARY, "sum_of_values");
		r[0].set(GRB_DoubleAttr_UB, 0);

		// Forbid edge from node back to itself
		for (int i = 0; i <= n; i++){
			x[i][i].set(GRB_DoubleAttr_UB, 0);
			x[i][0].set(GRB_DoubleAttr_UB, 0);
		}

		for (int v = 1; v <= n; v++){
			expr = 0;
			for (auto u : G.adj[v])
				expr += x[v][u.first];
			model.addConstr(expr >= y[v], "in_tree_terminal"+itos(v));
		}

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
			double *sol = model.get(GRB_DoubleAttr_X, y, n + 1);
			
			double **edge = new double*[n + 1];
			for (int v = 1; v <= n; v++)
				edge[v] = model.get(GRB_DoubleAttr_X, x[v], n + 1);

			cout << "Tree: ";
			for (int i = 1; i <= n; i++)
				cout << sol[i] << " ";
			cout << endl;
			
			for (int v = 1; v <= n; v++)
				for (int u = 1; u < v; u++){
					if (edge[u][v] > 0)
						cout<<u<<' '<<v<<' '<<edge[u][v]<<endl;
				}

			delete[] edge;			
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
