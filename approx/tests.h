#ifndef TESTS
#define TESTS
#include <sstream>

string itos(int i) {stringstream s; s << i; return s.str(); }

map<string, int> opt_sol;

void read_sol(){
	ifstream file;
	file.open("../optimal_solutions");

	if (!file.is_open())
		return;

	opt_sol.clear();

	int val;
	string s;
	while (file>>s>>val)
		opt_sol[s] = val;
}

int ninst, nopt;
double sumgap, mingap, maxgap;

void int_statistics(){
	ninst = 0;
	nopt = 0;
	sumgap = 0;
	mingap = 100;
	maxgap = 0;
}

void print_statistics(){
	cout<<"min "<<mingap<<" max "<<maxgap<<" avg "<<(sumgap/ninst)<<" nopt "<<nopt<<endl;
}


void execute(Graph G, string nme){
	double val = pcst(G).val();
	double opt = opt_sol[nme];
	double frac = ((val - opt)/val * 100.0);
	cout<<"nome    valor   otimo     gap%"<<endl;
	cout<<nme<<' '<<val<<' '<<opt<<' '<<frac<<"%"<<endl;
	mingap = min(mingap, frac);
	maxgap = max(maxgap, frac);
	sumgap += frac;
	ninst++;
	if (fabs(opt_sol[nme] - val) < 1) nopt++;
}



void runp(){
	int_statistics();
	for (int i = 0; i < 11; i++){
		string nme = "P100";
		if (i > 0)
			nme += "." + itos(i);
			
		if (i > 4)
			nme = "P200";
		if (i > 5){
			nme = "P400";
			if (i > 6)
				nme += "." + itos(i - 6);
		}

		string file_name = "../instances/PCSPG-JMP/" + nme  +".stp";
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}
		execute(G, nme);
	}
	print_statistics();
}	

void runk(){
	int_statistics();
	for (int i = 0; i < 23; i++){
		string nme = "K100";
		if (i > 0)
			nme += "." + itos(i);
		
		if (i > 10)
			nme = "K200";
		if (i > 11){
			nme = "K400";
			if (i > 12)
				nme += "." + itos(i - 12);
		}

		string file_name = "../instances/PCSPG-JMP/" + nme  +".stp";
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	print_statistics();
}	

void runc(){
	int_statistics();
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
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	print_statistics();
}	

void rund(){
	int_statistics();
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
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	print_statistics();
}	

void rune(){
	int_statistics();
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
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	print_statistics();
}

void runall(){
	runp();
	runk();
	runc();
	rund();
	rune();
}

#endif