/* Author : Gabriel Morete de Azevedo
   Feofiloff et al implementation of the JMP algorithm
   Union-Find library
*/

#ifndef DSU
#define DSU

#include <vector>
using namespace std;
struct Dsu {
	vector<int> id, sz;

	Dsu(int n) : id(n + 1), sz(n + 1) {
		for (int itr = 0; itr < n + 1; itr++){
			id[itr] = itr;
			sz[itr] = 1;
		} 
	};

	int find(int a){
		if (id[a] == a) return a;
		return id[a] = find(id[a]);
	}

	bool merge(int a, int b){
		a = find(a);
		b = find(b);
		if (a == b)
			return false;
		if (sz[b] > sz[a])
			swap(a, b);
		id[b] = a;
		sz[a] += sz[b];
		return true;
	}
};

#endif