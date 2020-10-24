/* Author : Gabriel Morete de Azevedo
   STP reader
*/

#include<iostream>
#include<fstream>
#include<vector>
using namespace std;

// Basic graph structure
struct Graph {
	int V, E;
	vector<int> p;
	vector< vector<pair<int, int> > > adj;
	Graph() {} Graph(int _V, int _E) : V(_V), E(_E), adj(V + 1), p(V + 1, 0) {};
};

// STP file format reader
// Returns the read Graph in G
// Returns 0 (sucess), 1 (error reading file), 2 (wrong format)
int STP_reader(string file_name, Graph &G){
	ifstream stp_file;

	stp_file.open(file_name);

	if (!stp_file.is_open())
		return 1;
	
	string s, coment;

	getline(stp_file, coment); 	// get comment
								// not a lot of uses for it now
	int V, E, cnt, a, b, c, nt;
	while (!stp_file.eof()){	
		stp_file >> s;
		if (s == "EOF") 
			break;
		if (s != "SECTION")
			return 2;
		
		stp_file >> s;
		if (s == "Comments")
			while (s != "END")
				stp_file>>s;	// do correct later, not much use now
		else if (s == "Graph"){
			for (int i = 0; i < 2; i++){
				stp_file >> s;
				if (s == "Nodes")
					stp_file >> V;
				if (s == "Edges")
					stp_file >> E;
			}

			G = Graph(V, E);

			stp_file>>s;
			while (s != "END"){
				stp_file >> a >> b >> c;
				G.adj[a].push_back({b, c});
				if (s == "E")
					G.adj[b].push_back({a, c});
				stp_file >> s;
				E--;
			}
			if (E)
				return 2;
		}	
		else if (s == "Terminals"){
			stp_file >> s >> nt;
			stp_file >> s;
			while (s != "END"){
				stp_file >> a >> c;
				G.p[a] = c;
				stp_file>>s;
				nt--;
			}
			if (nt)
				return 2;
		}
	}

	stp_file.close();
	return 0;
}

// int main(){
// 	string s = "K100.1.stp";
// 	Graph G;
// 	dbg(STP_reader(s, G));
// }