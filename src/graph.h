#ifndef GRAPH
#define GRAPH

#include<vector>
using namespace std;
// Basic graph structure
struct Graph {
	int V, E;
	vector<int> p;
	vector< vector<pair<int, int> > > adj;
	Graph() {} Graph(int _V, int _E) : V(_V), E(_E), adj(V + 1), p(V + 1, 0) {};
};

#endif 