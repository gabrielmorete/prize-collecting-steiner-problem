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

#define PCST // Signal for auxiliary libraries

#include "gurobi_c++.h"
#include "stp_reader.h"
#include "preprocessing.cpp"
#include "debug.h"
#include "dinic.h"
#include <cassert>
#include <cstdlib>	
#include <cstring>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

class cut_tree: public GRBCallback{
	public:
		GRBVar* vary; 
		GRBVar** varx;
		vector<vector<bool>> adj;
		int n;
		bool *term;
		cut_tree(GRBVar* _y, GRBVar** _x, int _n, vector<vector<bool>> _adj, bool * _term) {
			vary =_y; varx = _x;  n = _n; adj = _adj; term = _term;
		}
	protected:
		void callback() {
			try {
				if (where == GRB_CB_MIPSOL) {
					// Found an integer feasible solution - does it form a tree?
					//double dd = getDoubleInfo(GRB_CB_MIP_OBJBST);
					//	cout<<getDoubleInfo(GRB_CB_MIP_OBJBST)<<endl;

					double y[n + 1];
					vector<vector<double>> x(n + 1, vector<double>(n + 1));
					
					for (int v = 0; v <= n; v++)
						y[v] = getSolution(vary[v]);

					for (int v = 0; v <= n; v++){
						for (int u = 0; u <= n; u++){
							x[v][u] = getSolution(varx[v][u]);
					//		if (x[v][u] > 0)
					//		dbg(x[v][u]);
						}
					}

					for (int vtx = 1; vtx <= n; vtx++){
						if (!term[vtx])
							continue;

						init(); // start dinic	

						for (int v = 0; v <= n; v++)
							for (int u = 0; u <= n; u++)
								if (x[v][u] > 0){
									add(v, u, x[v][u]);
								}
		
						int cnt = 1;
						double f = dinic(0, vtx);		
						while (sign(f - y[vtx]) < 0 and cnt--){
							
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

	
							for (int v = 0; v <= n; v++){
								if (dist[v] >= 0){ // terminal inside the cut
									for (int u = 1; u <= n; u++){
										if (dist[u] < 0 and adj[v][u])
											add(v, u, 1);
									}
								}	
							}

							f += dinic(0, vtx);
						}
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
	GRBEnv env = GRBEnv(); // Gurobi env
	
	for (int i = 13; i < 41; i++){
		//string nme = "P400." + itos(i);
		string nme = "D" + itos(i) + "-A";
		if (i < 10)
			nme = "D0" + itos(i) + "-A";
		if (i > 20){
			nme = "D" + itos(i - 20) + "-B";
			if (i < 30)
				nme = "D0" + itos(i - 20) + "-B";
		}

		// if (i == 0)
		// 	nme = "P400";
//		string file_name = "../instances/PCSPG-JMP/" + nme  +".stp";
		string file_name = "../instances/PCSPG-CRR/" + nme  +".stp";
		//dbg(file_name);
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}

		int mold = G.E;

		// Reduction steps
		int en = preprocessing(G);
	
		double pn = 100*((double) (G.V - en)/G.V);
		double pm = 100*((double) (mold - G.E)/mold);

		int n = G.V, m = G.E;
		bool term[n + 1];
		vector<vector<bool>> adj(n + 1, vector<bool>(n + 1, 0));
		dtype cost[n + 1][n + 1];

		// printf("%s %d/%d(%.2f) %d/%d(%.2f)\n", nme.c_str(), n, en, pn, mold, m, pm);
		// fflush(stdout);	
		// continue;


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

		GRBVar **x = new GRBVar*[n + 1];
		for (int i = 0; i <= n; i++)
			x[i] = new GRBVar[n + 1];

		GRBVar *y = new GRBVar[n + 1];
		GRBVar *r = new GRBVar[n + 1];

		try {

			GRBModel model = GRBModel(env);

			// We will use an event Handler
			// Setting some parameters
			model.set(GRB_IntParam_OutputFlag, 0);
			model.set(GRB_IntParam_LazyConstraints, 1);
			model.set(GRB_DoubleParam_TimeLimit, 2000.0);
			model.set(GRB_IntParam_Presolve, 2);
		//	model.set(GRB_DoubleParam_Heuristics, 0.15);
			

			// Create binary decision variables
			for (int v = 0; v <= n; v++)
				for (int u = 0; u <= n; u++)
					x[v][u] = model.addVar(0.0, (double) adj[v][u], 0, GRB_CONTINUOUS, "x_"+itos(v)+"_"+itos(u));
			
			y[0] = model.addVar(0.0, 1.0, 0, GRB_BINARY, "y_"+itos(0));
			for (int v = 1; v <= n; v++)
				y[v] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "y_"+itos(v));
			
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

			// Exclude inexisting edges
			for (int v = 1; v <= n; v++)
				for (int u = 1; u <= v; u++)
					if (!adj[v][u])
						x[v][u].set(GRB_DoubleAttr_UB, 0);

			// Exclude inexisting verticies
			for (int v = 1; v <= n; v++)
				if (!term[v] and adj[v].size() == 0)		
					y[v].set(GRB_DoubleAttr_UB, 0);

			// Set Obj
			dtype obj = 0;
			for (int v = 1; v <= n; v++)
				obj += G.p[v];
			expr = 0;
			for (int v = 0; v <= n; v++)
				for (int u = 0; u <= n; u++)
					if (adj[v][u]) 
						expr += cost[v][u] * x[v][u];

			expr +=	obj;

			model.setObjective(expr, GRB_MINIMIZE);	

			// Set callback function
			cut_tree cb = cut_tree(y, x, n, adj, term);
			model.setCallback(&cb);

			// Optimize model
			model.optimize();

			// Extract solution

			if (model.get(GRB_IntAttr_SolCount) > 0) { /////////////////////////
				double *sol = model.get(GRB_DoubleAttr_X, y, n + 1);
				
				double **edge = new double*[n + 1];
				for (int v = 0; v <= n; v++)
					edge[v] = model.get(GRB_DoubleAttr_X, x[v], n + 1);

				// cout << "Tree: ";
				// for (int i = 0; i <= n; i++)
				// 	if(sol[i] >= 0.5)
				// 		cout<<i<<' ';
				// cout << endl;
				// cout<<"Edges"<<endl;
				// for (int v = 0; v <= n; v++){
				// 	for (int u = 0; u <= n; u++){
				// 		 if (edge[u][v] > 0.5)
				// 			cout<<u<<' '<<v<<' '<<cost[u][v]<<endl;
				// 	}
				// }	

				double tme = model.get(GRB_DoubleAttr_Runtime);
				dtype opt = (dtype) model.get(GRB_DoubleAttr_ObjVal);
				printf("%s %d/%d(%.2f) %d/%d(%.2f) %.3f %d\n", nme.c_str(), n, en, pn, mold, m, pm, tme, opt);
				fflush(stdout);	
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
	}
	return 0;
}

// void callback() {
// 	try {
// 		if (where == GRB_CB_MIPSOL) {
// 			// Found an integer feasible solution - does it form a tree?
// 			//cout<<getDoubleInfo(GRB_CB_MIP_OBJBST)<<endl;

// 			double y[n + 1];
// 			vector<vector<double>> x(n + 1, vector<double>(n + 1));
			
// 			for (int v = 0; v <= n; v++)
// 				y[v] = getSolution(vary[v]);

// 			for (int v = 0; v <= n; v++){
// 				for (int u = 0; u <= n; u++)
// 					x[v][u] = getSolution(varx[v][u]);
// 			}

// 			init(); // start dinic	

// 			for (int v = 0; v <= n; v++)
// 				for (int u = 0; u <= n; u++)
// 					if (x[v][u] > 0){
// 						add(v, u, x[v][u]);
// 					}

// 			for (int vtx = 1; vtx <= n; vtx++){
// 				if (sign(dinic(0, vtx) - y[vtx]) < 0){
// 					GRBLinExpr cut = 0;

// 					for (int v = 0; v <= n; v++){
// 						if (dist[v] >= 0){ // terminal inside the cut
// 							for (int u = 1; u <= n; u++){
// 								if (dist[u] < 0)
// 									cut += varx[v][u];
// 							}
// 						}	
// 					}

// 					addLazy(cut >= vary[vtx]);
// 				}
		
// 			}
// 		}
// 		} catch (GRBException e) {
// 			cout << "Error number: " << e.getErrorCode() << endl;
// 			cout << e.getMessage() << endl;
// 		} catch (...) {
// 			cout << "Error during callback" << endl;
// 		}
// }		