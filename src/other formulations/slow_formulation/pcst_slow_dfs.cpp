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

	Here we will try to use a simpler event handler (a dfs from the root) to
	find infeaseble cuts.

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
					// Found an integer feasible solution - does it form a tree?

					double y[n + 1], r[n + 1];
					vector<vector<double>> x(n + 1, vector<double>(n + 1));
					
					for (int i = 1; i <= n; i++){
						y[i] = getSolution(vary[i - 1]);
						r[i] = getSolution(varr[i - 1]);
					}

					for (int i = 1; i <= n; i++){
						for (int j = 1; j <= n; j++)
							x[i][j] = getSolution(varx[i - 1][j - 1]);
					}

					int root = 1, cnt = 0;
					for (int v = 2; v <= n; v++)
						if (r[v] > r[root])
							root = v;

					int tree_nodes;	

					// cout<<"Callback "<<callbackcnt++<<endl;	
					// cout<<"vtx : ";
					for (int v = 1; v <= n; v++)
						if (y[v] > 0.5){
//							cout<<v<<' ';
							tree_nodes++;
						}
					// gnl;
					// cout<<"root : "<<root<<endl;			

					int e = 0;
					for (int v = 1; v <= n; v++){
						for (int u = 1; u < v; u++)
							if (x[u][v] > 0){
//								cout<<u<<' '<<v<<' '<<x[u][v]<<endl;
								e++;
							}
					}

					// cout<<"#edges : "<<e<<endl;


					int vis[n + 1];
					memset(vis, 0, sizeof vis);

					if (dfs(root, n, vis, y, x) != tree_nodes){
			//			cout<<"Adicionado corte"<<endl;
						
						GRBLinExpr lhs = 0, rhs = 0;
						
			//			cout<<"Cut : ";		
						double cut_size = 0;		
						for (int v = 1; v <= n; v++)
							if (vis[v] == 0){
			//					cout<<v<<' ';
								rhs += vary[v - 1];
								cut_size++;
							}
							
				//		cout<<"("<<cut_size<<")"<<endl;


						for (int v = 1; v <= n; v++)
							if (vis[v] == 0){
				//				cout<<v<<' '<<r[v]<<endl;
								if (r[v] > 0.5){
									cout<<"Raiz no corte!"<<endl;
									assert(0);
								}

								rhs -= cut_size * varr[v - 1];
							}

						lhs = 0;
				//		cout<<"begin edges"<<endl;
					
						for (int v = 1; v <= n; v++){
							if (vis[v] == 0){ // sou terminal no corte
								for (int u = 1; u <= n; u++){
									if (vis[u] == 1){
										lhs += cut_size * varx[v - 1][u - 1];
				//						cout<<v<<' '<<u<<endl;
									}
								}
							}	
						}
				//		cout<<"end edges"<<endl;

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

	string file_name = "ljubic1.stp";
//	string file_name = "B/b03.stp";


	int _code = STP_reader(file_name, G);

	if (_code != 0){
		cout<<"Error reading file - Code "<<_code<<endl;
		return 0;
	}

	int n = G.V, m = G.E;
	double adj[n][n];

	int inf = 0;

	for (int v = 1; v <= n; v++){ // grafo esta indexado por um
		for (auto u : G.adj[v])  // variaveis por zero
			inf += adj[v - 1][u.first - 1];
		inf += G.p[v];
	}
	inf++;	

	for (int v = 0; v < n; v++)
		for (int u = 0; u < n; u++)
			adj[v][u] = inf;

	for (int v = 0; v < n; v++)	
		adj[v][v] = 0;

	for (int v = 1; v <= n; v++) // grafo esta indexado por um
		for (auto u : G.adj[v])  // variaveis por zero
			adj[v - 1][u.first - 1] = u.second;

	GRBEnv env = GRBEnv();
	GRBVar **x = new GRBVar*[n];
	for (int i = 0; i < n; i++)
		x[i] = new GRBVar[n];

	GRBVar *y = new GRBVar[n];
	GRBVar *r = new GRBVar[n];

	try {

		GRBModel model = GRBModel(env);

		// We will use an event Handler

		model.set(GRB_IntParam_LazyConstraints, 1);
		model.set(GRB_DoubleParam_TimeLimit, 1000.0);

		// Create binary decision variables

		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
				x[j][i] = x[i][j];
			}
		}

		GRBLinExpr expr = 0;
		for (int v = 0; v < n; v++){
			y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_"+itos(v));
			r[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_"+itos(v));
		}

		// Forbid edge from node back to itself
		for (int i = 0; i < n; i++)
			x[i][i].set(GRB_DoubleAttr_UB, 0);


		// Root is in the tree
		for (int v = 0; v < n; v++)
			model.addConstr(y[v] >= r[v], "in_tree_root_"+itos(v));

		// Only one root
		expr = 0;
		for (int v = 0; v < n; v++)	
			expr += r[v];
		model.addConstr(expr == 1, "one_root");

		// Incident inside the tree
		for (int v = 0; v < n; v++){
			expr = 0;
			for (int u = 0; u < n; u++)
				expr += x[v][u];
			model.addConstr(n * y[v] >= expr, "in_edge_in_tree"+itos(v));
		}

		// Set Obj
		expr = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < i; j++) 
				expr += adj[i][j] * x[i][j];

		for (int v = 0; v < n; v++)
				expr +=	G.p[v + 1] * ( 1 - y[v]);

		model.setObjective(expr, GRB_MINIMIZE);	

		// Set callback function

		cut_tree cb = cut_tree(y, r, x, n);
		model.setCallback(&cb);

		// Optimize model

		model.optimize();

		// Extract solution

		if (model.get(GRB_IntAttr_SolCount) > 0) { /////////////////////////
			double *sol = model.get(GRB_DoubleAttr_X, y, n);
			
			double **edge = new double*[n];
			for (int v = 0; v < n; v++)
				edge[v] = model.get(GRB_DoubleAttr_X, x[v], n);

			cout << "Tree: ";
			for (int i = 0; i < n; i++)
				cout << (sol[i] < 0.5 ? 0 : 1) << " ";
			cout << endl;
			
			for (int v = 0; v < n; v++){
				for (int u = 0; u < v; u++){
					 if (edge[u][v] > 0.5)
						cout<<u + 1<<' '<<v + 1<<' '<<adj[u][v]<<endl;
				}
			}	

			for (int v = 0; v < n; v++)
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

	// fr(i, n){
	// 	fr(j, n)
	// 		cout<<adj[i][j]<<' ';
	// 	gnl;		
	// }


	for (int i = 0; i < n; i++)
		delete[] x[i];
	delete[] x;
	delete[] y;
	delete[] r;

	return 0;
}
