/* Author : Gabriel Morete de Azevedo
   Branch and Bound algorithm
*/

#include<iostream>
#include<fstream>
#include<vector>
#include"debug.h"
using namespace std;

struct node{
	int lb, ub, n_edge, n_vtx;
	vector< pair<int, int> > edge, vertex;	
	node() {} node(int _lb, int _ub) : lb(_lb), ub(_ub), edge(), p() {};

};

int branch_and_Bound(Graph G){
	int lb = 0, ub = inf;
	vector<nodes> active_list;
	active_list.pb(node(lb, ub));

	while (!active_list.empty()){
		
	}

}