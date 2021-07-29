/* Author : Gabriel Morete de Azevedo
   Preprocessing library for PCST
*/

#include "prepro.h"
#include "stp_reader.h"
#include "graph.h"
#include "pcdist.cpp"
#include "dsu.h"
#include "debug.h"
#include <map>
#include <chrono>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

typedef vector< vector< pair< int, type_val> > > vvit;

bool mst_test(int &n, int &m, vector<bool> &term, vvit &adj, vector<type_val> &p){
	bool tried[n + 1]; // tentei o vértice	
	int pnew[n + 1]; // novo p
	int min_cst[n + 1]; // menor custo
	Dsu dsu = Dsu(n); // union find

	vector< tuple<int, int, int> > edges;
	for (int v = 1; v <= n; v++)
		for (auto e : adj[v])
			if (v < e.first)
				edges.push_back({e.second, v, e.first}); 

	sort(edges.begin(), edges.end()); // arestas ordenadas por custo	

	for (int v = 1; v <= n; v++)
		pnew[v] = p[v];

	bool ok = 1;
	while (ok){
		ok = 0;
		memset(tried, 0, sizeof tried);

		for (auto e : edges){
			int v, u, cs;
			tie(cs, v, u) = e;
			int rv = dsu.find(v); // representante de v
			int ru = dsu.find(u); // representante de u
			if (min(pnew[rv], pnew[ru]) > cs and ru != rv){ 

				if ((tried[rv] and cs != min_cst[rv]) or (tried[ru] and cs != min_cst[ru]))
					continue; // aresta não tem custo minimo

				dsu.merge(u, v);
				
				pnew[dsu.find(v)] = pnew[rv] + pnew[ru] - cs;
				ok = 1;
				break;
			}
			else if (ru != rv){ // não adicionei aresta, agora tem cuso minimo
				if (!tried[rv]){
					tried[rv] = 1;
					min_cst[rv] = cs;
				}
				if (!tried[ru]){
					tried[ru] = 1;
					min_cst[ru] = cs;
				}
			}
		}
	}	

	vector< tuple<int, int, int> > edges2;

	for (auto e : edges) // aresta cruza componentes
		if (dsu.find(get<1>(e)) != dsu.find(get<2>(e)))
			edges2.push_back(e);

	int removed = 0;	

	map< pair<int, int>, int > join; // não posso usar dsu

	for (auto e : edges2){
		int v, u, cs;
		tie(cs, v, u) = e;	

		if (!join.count({dsu.find(v), dsu.find(u)})){
			join[{dsu.find(v), dsu.find(u)}] = 1;
			join[{dsu.find(u), dsu.find(v)}] = 1;
		}
		else {
			m--;
			removed++;
		}
	}	

	// if(removed != 0)
	// 	cout<<"mst : "<<removed<<endl;

	return removed > 0;
}

void print(string s, int mo, int m){
	cout<<s<<": "<<mo<<' '<<m<<' '<<100*((double) (mo - m)/mo)<<"%"<<endl;
}


// Preprocessing function, if some error occour the return value is -1,
// othwewise is the number of real vertices remaining in the graph.
void preprocessing(Graph &G){
	int n = G.V;
	int m = G.E;
	vector<type_val> p = G.p;
	vvit adj = G.adj;

	vector<bool> term(n + 1);
	for (int v = 1; v <= n; v++)
		term[v] = G.p[v] > 0;


	mst_test(n, m, term, adj, p);	
	
	G.adj = adj;
	G.E = m;
}

void run(vector<string> names, string folder){
	Graph G;
	for (auto nme : names){
		string file_name = "../instances/PCSPG-" + folder + "/" + nme  +".stp";
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);
	}
}


#ifndef PCST

int main(){
	// string path = "../instances/";

	// string s = "E/E20-A.stp";

	Graph G;
	
	// cout<<setprecision(2)<<fixed;

	// if (STP_reader(path + s, G) != 0){
	// 	cout<<"Arquivo não existe"<<endl;
	// 	return 0;
	// }

	// int mo = G.E;
	// preprocessing(G);
	// int m = G.E;	




	for (int i = 0; i < 10; i++){
		string nme = "P100";
		if (i > 0)
			nme += "." + itos(i);
			

		if (i > 4)
			nme = "P200";
		if (i > 5){
			nme = "P400";
			if (i > 6)
				nme += "." +itos(i - 5);
		}

		string file_name = "../instances/PCSPG-JMP/" + nme  +".stp";
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}

		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);

	}

	for (int i = 0; i < 23; i++){
		string nme = "K100";
		if (i > 0)
			nme += "." + itos(i);
		
		if (i > 10)
			nme = "K200";
		if (i > 11){
			nme = "K400";
			if (i > 12)
				nme += "." +itos(i - 12);
		}

		string file_name = "../instances/PCSPG-JMP/" + nme  +".stp";
		
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}
		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);

	}

	for (int i = 1; i < 41; i++){
		string nme = "C" + itos(i) + "-A";
		if (i < 10)
			nme = "C0" + itos(i) + "-A";
		if (i > 20){
			nme = "C" + itos(i - 20) + "-B";
			if (i < 30)
				nme = "C0" + itos(i - 20) + "-B";
		}

		string file_name = "../instances/PCSPG-CRR/" + nme  +".stp";
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}
		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);
	}


	for (int i = 1; i < 41; i++){
		string nme = "D" + itos(i) + "-A";
		if (i < 10)
			nme = "D0" + itos(i) + "-A";
		if (i > 20){
			nme = "D" + itos(i - 20) + "-B";
			if (i < 30)
				nme = "D0" + itos(i - 20) + "-B";
		}

		string file_name = "../instances/PCSPG-CRR/" + nme  +".stp";
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}
		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);
	}

	for (int i = 1; i < 41; i++){
		string nme = "E" + itos(i) + "-A";
		if (i < 10)
			nme = "E0" + itos(i) + "-A";
		if (i > 20){
			nme = "E" + itos(i - 20) + "-B";
			if (i < 30)
				nme = "E0" + itos(i - 20) + "-B";
		}

		string file_name = "../instances/E/" + nme  +".stp";
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return 0;
		}
		int mo = G.E;
		preprocessing(G);
		int m = G.E;	
		print(nme, mo, m);
	}

	vector<string> H = {"hc10p", "hc11p", "hc12p", "hc6p", "hc7p", "hc8p", "hc9p",
	"hc10u", "hc11u", "hc12u", "hc6u", "hc7u", "hc8u", "hc9u"};
	// run(H, "H");

	vector<string> PUCNU = {
	"bip42nu",  "bipa2nu" ,  "cc11-2nu",  "cc3-11nu",  "cc3-5nu",  "cc6-3nu",
	"bip52nu",  "bipe2nu" ,  "cc12-2nu",  "cc3-12nu",  "cc5-3nu",  "cc7-3nu",
	"bip62nu",  "cc10-2nu",  "cc3-10nu",  "cc3-4nu" ,  "cc6-2nu",  "cc9-2nu"};
	// run(PUCNU, "PUCNU");

	vector<string> RAND = {
	"a0200RandGraph.1.2",   "a1000RandGraph.3",     "a1800RandGraph.2",
	"a0200RandGraph.1.5",   "a12000RandGraph.1.2",  "a1800RandGraph.3",
	"a0200RandGraph.2",     "a12000RandGraph.1.5",  "a2000RandGraph.1.2",
	"a0200RandGraph.3",     "a12000RandGraph.2",    "a2000RandGraph.1.5",
	"a0400RandGraph.1.2",   "a12000RandGraph.3",    "a2000RandGraph.2",
	"a0400RandGraph.1.5",   "a1200RandGraph.1.2",   "a2000RandGraph.3",
	"a0400RandGraph.2",     "a1200RandGraph.1.5",   "a3000RandGraph.1.2",
	"a0400RandGraph.3",     "a1200RandGraph.2",     "a3000RandGraph.1.5",
	"a0600RandGraph.1.2",   "a1200RandGraph.3",     "a3000RandGraph.2",
	"a0600RandGraph.1.5",   "a14000RandGraph.1.2",  "a3000RandGraph.3",
	"a0600RandGraph.2",     "a14000RandGraph.1.5",  "a4000RandGraph.1.2",
	"a0600RandGraph.3",     "a14000RandGraph.2",    "a4000RandGraph.1.5",
	"a0800RandGraph.1.2",   "a14000RandGraph.3",    "a4000RandGraph.2",
	"a0800RandGraph.1.5",   "a1400RandGraph.1.2",   "a4000RandGraph.3",
	"a0800RandGraph.2",     "a1400RandGraph.1.5",   "a6000RandGraph.1.2",
	"a0800RandGraph.3",     "a1400RandGraph.2",     "a6000RandGraph.1.5",
	"a10000RandGraph.1.2",  "a1400RandGraph.3",     "a6000RandGraph.2",
	"a10000RandGraph.1.5",  "a1600RandGraph.1.2",   "a6000RandGraph.3",
	"a10000RandGraph.2",    "a1600RandGraph.1.5",   "a8000RandGraph.1.2",
	"a10000RandGraph.3",    "a1600RandGraph.2",     "a8000RandGraph.1.5",
	"a1000RandGraph.1.2",   "a1600RandGraph.3",     "a8000RandGraph.2",
	"a1000RandGraph.1.5",   "a1800RandGraph.1.2",   "a8000RandGraph.3",
	"a1000RandGraph.2",     "a1800RandGraph.1.5"};
	// run(RAND, "RANDOM");


	vector<string> ACTMOD = {
	"drosophila001",  "drosophila0075",  "lymphoma",             "metabol_expr_mice_2",
	"drosophila005",  "HCMV",            "metabol_expr_mice_1",  "metabol_expr_mice_3"};
	// run(ACTMOD, "ACTMODPC");

	vector<string> I640 ={
	"i640-001",  "i640-033",  "i640-115",  "i640-202",  "i640-234",  "i640-321",
	"i640-002",  "i640-034",  "i640-121",  "i640-203",  "i640-235",  "i640-322",
	"i640-003",  "i640-035",  "i640-122",  "i640-204",  "i640-241",  "i640-323",
	"i640-004",  "i640-041",  "i640-123",  "i640-205",  "i640-242",  "i640-324",
	"i640-005",  "i640-042",  "i640-124",  "i640-211",  "i640-243",  "i640-325",
	"i640-011",  "i640-043",  "i640-125",  "i640-212",  "i640-244",  "i640-331",
	"i640-012",  "i640-044",  "i640-131",  "i640-213",  "i640-245",  "i640-332",
	"i640-013",  "i640-045",  "i640-132",  "i640-214",  "i640-301",  "i640-333",
	"i640-014",  "i640-101",  "i640-133",  "i640-215",  "i640-302",  "i640-334",
	"i640-015",  "i640-102",  "i640-134",  "i640-221",  "i640-303",  "i640-335",
	"i640-021",  "i640-103",  "i640-135",  "i640-222",  "i640-304",  "i640-341",
	"i640-022",  "i640-104",  "i640-141",  "i640-223",  "i640-305",  "i640-342",
	"i640-023",  "i640-105",  "i640-142",  "i640-224",  "i640-311",  "i640-343",
	"i640-024",  "i640-111",  "i640-143",  "i640-225",  "i640-312",  "i640-344",
	"i640-025",  "i640-112",  "i640-144",  "i640-231",  "i640-313",  "i640-345",
	"i640-031",  "i640-113",  "i640-145",  "i640-232",  "i640-314",
	"i640-032",  "i640-114",  "i640-201",  "i640-233",  "i640-315"};
	// run(I640, "i640");

	vector<string> HBD = {
		"handbd01",  "handbd04",  "handbd07",  "handbd10",  "handbd13",
		"handbd02",  "handbd05",  "handbd08",  "handbd11",  "handbd14",
		"handbd03",  "handbd06",  "handbd09",  "handbd12"};
	// run(HBD, "hand/HAND_BIG_DIMACS");
	
	vector<string> HBI = {
		"handbi01",  "handbi04",  "handbi07",  "handbi10",  "handbi13",
		"handbi02",  "handbi05",  "handbi08",  "handbi11",  "handbi14",
		"handbi03",  "handbi06",  "handbi09",  "handbi12"};
	// run(HBI, "hand/HAND_BIG_ICERM");

	vector<string> HSD = {
		"handsd01",  "handsd03",  "handsd05",  "handsd07",  "handsd09",
		"handsd02",  "handsd04",  "handsd06",  "handsd08",  "handsd10"};
	// run(HSD, "hand/HAND_SMALL_DIMACS");

	vector<string> HSI = {
		"handsi01",  "handsi03",  "handsi05",  "handsi07",  "handsi09",
		"handsi02",  "handsi04",  "handsi06",  "handsi08",  "handsi10"};
	// run(HSI, "hand/HAND_SMALL_ICERM");
	

}

#endif
