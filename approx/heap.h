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

frac get_min_key(int typ, int p){
	int q = (*H[typ][p].begin()).first;
	return key(q, p);
}

int get_min_ele(int typ, int p){
	int q = (*H[typ][p].begin()).first;
	return q;
}

bool is_in(int typ, int p, int q){
	return H[typ][p].count({q, p}) > 0;
}

#endif