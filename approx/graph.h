/* Author : Gabriel Morete de Azevedo
   Basic graph struct
*/

#ifndef GRAPH
#define GRAPH

#include <vector>
using namespace std;
// Basic graph structure

using type_val = long long;

struct Graph {
	int V, E;
	vector<type_val> p;
	vector< vector<pair<int, type_val> > > adj;
	Graph() {} Graph(int _V, int _E) : V(_V), E(_E), adj(V + 1), p(V + 1, 0) {};
};

#endif 