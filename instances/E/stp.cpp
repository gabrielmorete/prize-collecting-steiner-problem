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
#define all(x)  x.begin(),x.end()

typedef long long int ll;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<pii> vii;

const int INF = 0x3f3f3f3f;
const ll llINF = (long long)(1e18) + 100;   
const int MAXN = 2500;

string itos(int i) {stringstream s; s << i; return s.str(); }

typedef tuple<int, int, int> tup;

void transform(string nme_in, string nme_out){
	std::ofstream out(nme_out);
    std::ifstream in(nme_in);

    std::cin.rdbuf(in.rdbuf());
    std::cout.rdbuf(out.rdbuf());   

	string s;

	getline(cin, s);
	getline(cin, s);

	// dbg(s);

	vector< pair<int, int> > terminal;

	int vtx, cst, a, b;
	for (int i = 1; i <= MAXN; i++){
		cin>>vtx>>a>>b>>cst;
		// cout<<vtx<<' '<<a<<' '<<b<<' '<<cst<<endl;
		if (cst > 0)
			terminal.push_back({vtx, cst});
	}

	getline(cin, s);
	getline(cin, s);
	getline(cin, s);

	// dbg(s);

	vector<tup> edges;

	int cnt;

	while (cin>>cnt>>a>>b>>cst)
		edges.push_back({a, b, cst});
	sort(all(edges));

	// Write file
	cout<<"STP File, STP Format Version 1.0"<<endl<<endl;
	
	cout<<"SECTION Comments"<<endl;
	cout<<"Name \""<<nme_out<<"\""<<endl;
	cout<<"Creator \"Lujibc and converted by Gabriel Morete\""<<endl;
	cout<<"Problem \"Prize-Collecting Steiner Problem in Graphs\""<<endl;
	cout<<"END"<<endl<<endl;

	cout<<"SECTION Graph"<<endl;
	cout<<"Nodes "<<MAXN<<endl;
	cout<<"Edges "<<edges.size()<<endl;
	for (auto e : edges){
		tie(a, b, cst) = e;
		cout<<"E "<<a<<' '<<b<<' '<<cst<<endl;
	}
	cout<<"END"<<endl<<endl;

	cout<<"SECTION Terminals"<<endl;
	cout<<"Terminals "<<terminal.size()<<endl;
	for (auto x : terminal)
		cout<<"TP "<<x.first<<' '<<x.second<<endl;

	cout<<"END"<<endl<<endl;

	cout<<"EOF"<<endl;
}

int32_t main(){
	for (int i = 1; i < 41; i++){
		string nme = "e" + itos(i) + ".stp-A";
		string nme2 = "E" + itos(i) + "-A.stp";

		if (i < 10){
			nme = "e0" + itos(i) + ".stp-A";
			nme2 = "E0" + itos(i) + "-A.stp";
		}
			
		if (i > 20){
			nme = "e" + itos(i - 20) + ".stp-B";
			nme2 = "E" + itos(i - 20) + "-B.stp";

			if (i < 30){
				nme = "e0" + itos(i - 20) + ".stp-B";
				nme2 = "E0" + itos(i - 20) + "-B.stp";
			}
		}

		string file_name = "original/" + nme;
	// 	ifstream stp_file;

	// 	stp_file.open(file_name);

	// 	dbg(file_name);
	// if (!stp_file.is_open())
	// 	cout<<"hhhhhhhhhh"<<endl;
		transform(file_name, nme2);
	}
}