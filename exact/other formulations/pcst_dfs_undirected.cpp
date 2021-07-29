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

#define PCST

#include "gurobi_c++.h"
#include "../stp_reader.h"
#include "../debug.h"
#include "../preprocessing.cpp"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.


const double EPS = 1e-8;
int singn(double x){ return (x > EPS) - (x < -EPS); }

typedef vector< vector<double> > matrix;	

// returns the # of visited nodes
int dfs(int v, int n, int *vis, double *y, matrix adj){
	vis[v] = 1;

	int tree = 1;
	for (int u = 1; u <= n; u++)
		if (y[u] > 0.5 and adj[v][u] > 0.5 and !vis[u])
			tree += dfs(u, n, vis, y, adj);

	return tree;	
}

int callbackcnt;

class cut_tree: public GRBCallback{
	public:
		GRBVar* vary; 
		GRBVar** varx;
		int n;
		cut_tree(GRBVar* _y, GRBVar** _x, int _n) {
			vary =_y; varx = _x;  n = _n;
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

					int tree_nodes;	

					for (int v = 0; v <= n; v++)
						if (y[v] > 0.5)
							tree_nodes++;

					int e = 0;
					for (int v = 0; v <= n; v++){
						for (int u = 0; u < v; u++)
							if (x[u][v] > 0)
								e++;
					}

					int vis[n + 1];
					memset(vis, 0, sizeof vis);

					if (dfs(0, n, vis, y, x) != tree_nodes){
						GRBLinExpr lhs = 0, rhs = 0;
						double cut_size = 0;		
						for (int v = 1; v <= n; v++)
							if (vis[v] == 0){
								rhs += vary[v];
								cut_size++;
							}

						lhs = 0;
					
						for (int v = 0; v <= n; v++){
							if (vis[v] == 1){ // sou terminal no corte
								for (int u = 1; u <= n; u++){
									if (vis[u] == 0){
										lhs += cut_size * varx[v][u];
									}
								}
							}	
						}

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

using dtype = long long int;

int main(){
	Graph G; // type graph from stp reader
	GRBEnv env = GRBEnv();

	for (int i = 0; i < 22; i++){
		//string nme = "ljubic" + itos(i);
		string nme = "K100." + itos(i);
		if (i == 0 or i == 11)
			nme = "K100";

		if (i > 11)
			nme = "K100." + itos(i - 11);


		string file_name = "../../instances/PCSPG-JMP/" + nme  +".stp";
		

		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}

		if (i > 10)
			preprocessing(G);

		int n = G.V, m = G.E;
		bool term[n + 1];
		dtype adj[n + 1][n + 1], dist[n + 1][n + 1];

		term[0] = 1;
		for (int v = 1; v <= n; v++){
			term[v] = 0;
			if (G.p[v] > 0)
				term[v] = 1;
		}

		int inf = -1;

		for (int v = 0; v <= n; v++)
			for (int u = 0; u <= n; u++)
				adj[v][u] = inf;

		for (int v = 0; v <= n; v++){	
			adj[v][v] = 0;
			adj[0][v] = 0;
		}	

		for (int v = 1; v <= n; v++) // graph is 1-indexed
			for (auto u : G.adj[v])  // root is node 0
				adj[v][u.first] = u.second;

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
			model.set(GRB_IntParam_OutputFlag, 0);
			model.set(GRB_IntParam_LazyConstraints, 1);
			model.set(GRB_DoubleParam_TimeLimit, 300.0);
			model.set(GRB_DoubleParam_Heuristics, 0.15);

			// Create binary decision variables
			for (int i = 0; i <= n; i++) {
				for (int j = 0; j <= i; j++) {
					x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
					x[j][i] = x[i][j];
				}
			}

			GRBLinExpr expr = 0;
			for (int v = 0; v <= n; v++)
				y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_"+itos(v));
			
			// Forbid edge from node back to itself
			for (int v = 0; v <= n; v++)
				x[v][v].set(GRB_DoubleAttr_UB, 0);

			// Add non-costumer verticies restiction
			for (int v = 1; v <= n; v++)
				if (!term[v]) // non costumer can't be root
					x[0][v].set(GRB_DoubleAttr_UB, 0);		

			// Root is in the tree
			y[0].set(GRB_DoubleAttr_LB, 1);

			// Only one root
			expr = 0;
			for (int v = 1; v <= n; v++)	
				expr += x[0][v];
			model.addConstr(expr == 1, "one_root");

			// Incident inside the tree
			for (int v = 1; v <= n; v++){
				expr = 0;
				for (int u = 0; u <= n; u++)
					expr += x[v][u];
				model.addConstr((n + 1) * y[v] >= expr, "in_edge_in_tree"+itos(v));
			}

			// Symmetry break
			for (int v = 1; v <= n; v++) // fix after non-costumer-restriction
				for (int u = 1; u < v; u++)
					if (term[u]) // non-costumer might mess
						model.addConstr(x[0][v] <= 1 - y[u], "symmetry_"+itos(v)+"_"+itos(u));

			// Tree size
			expr = 0;
			for (int v = 1; v <= n; v++)
				expr += y[v];
			for (int v = 1; v <= n; v++)
				for (int u = 1; u < v; u++)
					expr -= x[v][u];
			model.addConstr(expr == 1, "tree_size");					

			// Non-leaf non-costumer (Lucena)
			for (int v = 1; v <= n; v++){
				expr = 0;
				for (int u = 1; u <= n; u++)
					expr += x[v][u];
				if (!term[v])				
					model.addConstr(expr >= 2 * y[v]);
				// else // Non-single vertex
				// 	model.addConstr(expr >= y[v]);
			}	

			// Exclude inexisting edges
			for (int v = 1; v <= n; v++)
				for (int u = 1; u <= v; u++)
					if (adj[v][u] == inf)
						x[v][u].set(GRB_DoubleAttr_UB, 0);

			// Set Obj
			expr = 0;
			for (int i = 1; i <= n; i++)
				for (int j = 1; j < i; j++) 
					expr += adj[i][j] * x[i][j];

			for (int v = 1; v <= n; v++)
					expr +=	G.p[v] * ( 1 - y[v]);

			model.setObjective(expr, GRB_MINIMIZE);	

			// Set callback function
			cut_tree cb = cut_tree(y, x, n);
			model.setCallback(&cb);

			// Optimize model

			model.optimize();

			// Extract solution

			cout<<setprecision(2)<<fixed;

			if (model.get(GRB_IntAttr_SolCount) > 0) {
				// double *sol = model.get(GRB_DoubleAttr_X, y, n + 1);
				
				// double **edge = new double*[n + 1];
				// for (int v = 0; v <= n; v++)
				// 	edge[v] = model.get(GRB_DoubleAttr_X, x[v], n + 1);

				// cout << "Tree: ";
				// for (int i = 0; i <= n; i++)
				// 	if(sol[i] >= 0.5)
				// 		cout<<i<<' ';
				// cout << endl;
				// cout<<"Edges"<<endl;
				// for (int v = 0; v <= n; v++){
				// 	for (int u = 0; u <= v; u++){
				// 		 if (edge[u][v] > 0.5)
				// 			cout<<u<<' '<<v<<' '<<adj[u][v]<<endl;
				// 	}
				// }	
				// cout<<"preprocessing : "<<pren<<" "<<prem<<endl;		
				double tme = model.get(GRB_DoubleAttr_Runtime);
				cout<<nme<<' '<<tme<<' '<<(int)model.get(GRB_DoubleAttr_ObjVal)<<' '<<model.get(GRB_DoubleAttr_MIPGap) * 100<<endl;


				// for (int v = 0; v <= n; v++)
				// 	delete[] edge[v];
				// delete[] edge;			
				// delete[] sol;
			}
			else{
				cout<<nme<<' '<<model.get(GRB_DoubleAttr_MIPGap) * 100<<endl;
			}

		} catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		} catch (...) {
			cout << "Error during optimization" << endl;
		}

		// fr(i, n){
		// 	fr(j, n)
		// 		cout<<adj[i][j]<<' ';
		// 	gnl;		
		// }


		for (int i = 0; i <= n; i++)
			delete[] x[i];
		delete[] x;
		delete[] y;
		delete[] r;
	}

	return 0;
}
