#include "bits/stdc++.h"
using namespace std;

#define pb push_back
#define fst first
#define snd second

#define fr(i,n)     for (int i = 0; i < n; i++)
#define frr(i,n)    for (int i = 1; i <= n; i++)

#define endl '\n'
#define gnl cout << endl
#define chapa cout << "oi meu chapa" << endl

#define dbg(x)  cout << #x << " = " << x << endl
#define all(x)  x.begin(), x.end()

typedef long long int ll;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<pii> vii;

const int INF = 0x3f3f3f3f;
const ll llINF = (long long)(1e18) + 100;   
const int MAXN = 1e5 + 10;


// Heap vai guardar o elemento e comparar usando a função key

int key(int a){
	return a;
}

struct Heap {
	vector<int> v;
	int cap, sz, id;

	// Constructor
	Heap(int _c, int _id){
		cap = _c;
		sz = 0;
		id = _id;
		v.resize(_c);
	}	

    int pre(int i) { return (i - 1) / 2; }

    int left(int i) { return (2 * i + 1); }
  
    int right(int i) { return (2 * i + 2); }

	int getmin() { return v[0]; }

	bool less(int a, int b){ // comprator
		return a < b;
	}

	void insertkey(int k){
		if (sz == cap){
			cout << "\nOverflow: Could not insertKey\n";
			return;
		}

		int pos = sz++;
		v[pos] = k;

		while (pos != 0 and less(v[pos], v[pre(pos)])){ ///////////////////////////////
			swap(v[pos], v[pre(pos)]);
			pos = pre(pos);
		}
	}

	int extractmin(){
		if (sz <= 0)
			return INT_MAX;
		if (sz == 1){
			sz--;
			return v[0];
		}

		int root = v[0];
		v[0] = v[sz - 1];
		sz--;
		heapify(0);

		return root;
	}

	void heapify(int pos){
		int l = left(pos);
		int r = right(pos);
		int small = pos;
		if (l < sz and less(v[l], v[pos]))
			small = l;
		if (r < sz and less(v[r], v[small]))
			small = r;
		if (small != pos){
			swap(v[pos], v[small]);
			heapify(small);
		}
	}

	void decreasekey(int pos, int new_val){
		v[pos] = new_val;
		while (pos != 0 and less(v[pos], v[pre(pos)])){
			swap(v[pos], v[pre(pos)]);
			pos = pre(pos);
		}
	}

	void deletekey(int pos){
		decreasekey(pos, INT_MIN);
		extractmin();
	}

	int find(int val){
		int pos = 0;
		while (pos < sz){
			if (val == v[pos])
				return pos;
			if (less(val, v[pos]))
				pos = left(pos);
			else
				pos = right(pos);
			dbg(pos);
		}

		return -1;
	}

};


int32_t main(){
	ios_base::sync_with_stdio(false); cin.tie(NULL); cout.tie(NULL);
	Heap h = Heap(11, 4);
	
	h.insertkey(3);
	cout << h.getmin()<<endl;
	h.insertkey(2);
	cout << h.getmin()<<endl;
	// h.deletekey(1);
	// cout << h.getmin()<<endl;

	fr(i, 4)
		dbg(h.v[i]);

	dbg((h.find(3)));

	h.insertkey(15);
	h.insertkey(5);
	h.insertkey(4);
	h.insertkey(45);

	cout << h.extractmin() << " ";
	cout << h.getmin() << " ";
	h.decreasekey(2, 1);
	cout << h.getmin();
	gnl;
}