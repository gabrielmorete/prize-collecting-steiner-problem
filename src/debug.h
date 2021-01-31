#ifndef DEBUG
#define DEBUG

#define gnl cout << endl
#define chapa cout << "oi meu chapa" << endl

#define dbg(x)  cout << #x << " = " << x << endl
#define all(x)  x.begin(),x.end()

#define fr(i,n)     for (int i = 0; i < n; i++)
#define frr(i,n)    for (int i = 1; i <= n; i++)

const int INF = 0x3f3f3f3f;

void debug_adj(int v, vector<vector<pair<int, int>>> &adj){
	cout<<"*** Adj "<<v<<" ***"<<endl;
	for (auto x : adj[v])
		cout<<x.first<<' '<<x.second<<endl;

	cout<<"*******"<<endl;
}


#endif 