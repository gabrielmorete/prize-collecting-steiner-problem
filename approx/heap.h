/* Author : Gabriel Morete de Azevedo
   Auxiliary functions to use std::set as a heap
*/

#ifndef HEAP
#define HEAP

bool insert(int typ, int p, int q){
	H[typ][p].insert({q, p});
	return true;
}

bool remove(int typ, int p, int q){
	if (H[typ][p].count({q, p})){
		H[typ][p].erase({q, p});
		return true;
	}
	return false;
}

void copy(int typ, int p, int q){
	for (auto x : H[typ][p])
		H[typ][q].insert({x.first, q});
}

int get_min_ele(int typ, int p){
	int q = (*H[typ][p].begin()).first;
	if (!A.count({p, q})){
		cout << "\033[1;31mRuntime error. Inconsistent Inner Structure. Aborting!\033[0m\n";
		exit(-1);
	}
	return q;
}

frac get_min_key(int typ, int p){
	int q = get_min_ele(typ, p);
	if (key(p, q) < frac(0)){
		cout << "\033[1;31mRuntime error. Violated Edge. Aborting!\033[0m\n";
		exit(-1);	
	}
	return key(p, q);
}

bool is_in(int typ, int p, int q){
	return H[typ][p].find({q, p}) != H[typ][p].end();
}

#endif