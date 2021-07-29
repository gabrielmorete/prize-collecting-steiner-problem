/* Author : Gabriel Morete de Azevedo
   Simple library that runs automated tests
*/

#ifndef TESTS
#define TESTS
#include <iomanip>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

map<string, int> opt_sol;

void read_sol(){ // read solution values
	ifstream file;
	file.open("../instances/optimal_values");

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
	cout<<setprecision(2)<<fixed;
	cout<<"min "<<mingap<<" max "<<maxgap<<" avg "<<(sumgap/ninst)<<" nopt "<<nopt<<endl;
}


void execute(Graph G, string nme){
	double val = pcst(G).val();
	double opt = opt_sol[nme];
	double frac = ((val - opt)/opt * 100.0);
	cout<<setprecision(2)<<fixed;
	cout<<nme<<' '<<val<<' '<<opt<<' '<<frac<<"%"<<endl;
	mingap = min(mingap, frac);
	maxgap = max(maxgap, frac);
	sumgap += frac;
	ninst++;
	if (fabs(opt_sol[nme] - val) < 1) nopt++;
}



void runp(){
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

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
	cout<<"P  ";
	print_statistics();
}	

void runk(){
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

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
	cout<<"K  ";
	print_statistics();
}	

void runc(){
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

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
	cout<<"C  ";
	print_statistics();
}	

void rund(){
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

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
	cout<<"D  ";
	print_statistics();
}	

void rune(){
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

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
	cout<<"E  ";
	print_statistics();
}

void runh(){
	vector<string> names = {"hc6p",	"hc6u", "hc7p", "hc7u", "hc8p", "hc8u",
		"hc9p", "hc9u", "hc10p", "hc10u", "hc11p", "hc11u", "hc12p" };
	
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

	for (string nme : names){

		string file_name = "../instances/PCSPG-H/" + nme + ".stp";
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	cout<<"H  ";
	print_statistics();
}

void runi640(){
	vector<string> names = {"i640-202",  "i640-234",  "i640-321",
		"i640-203",  "i640-235",  "i640-322",
		"i640-204",  "i640-241",  "i640-323",
		"i640-205",  "i640-242",  "i640-324",
		"i640-211",  "i640-243",  "i640-325",
		"i640-212",  "i640-244",  "i640-331",
		"i640-213",  "i640-245",  "i640-332",
		"i640-214",  "i640-301",  "i640-333",
		"i640-215",  "i640-302",  "i640-334",
		"i640-221",  "i640-303",  "i640-335",
		"i640-222",  "i640-304",  "i640-341",
		"i640-223",  "i640-305",  "i640-342",
		"i640-224",  "i640-311",  "i640-343",
		"i640-225",  "i640-312",  "i640-344",
		"i640-231",  "i640-313",  "i640-345",
		"i640-232",  "i640-314",  "i640-201",
		"i640-233",  "i640-315"};

	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

	for (string nme : names){

		string file_name = "../instances/PCSPG-i640/" + nme + ".stp";
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	cout<<"i640  ";
	print_statistics();
}

void runpcnu(){
	vector<string> names = {
		"bip42nu", "bip52nu", "bip62nu", "bipa2nu",
		"bipe2nu", "cc10-2nu", "cc11-2nu", "cc12-2nu",
		"cc3-10nu", "cc3-11nu", "cc3-12nu", "cc3-4nu",
		"cc3-5nu", "cc5-3nu", "cc6-2nu", "cc6-3nu",
		"cc7-3nu","cc9-2nu" };
	int_statistics();
	cout<<"nome    valor    otimo    gap(%)"<<endl;

	for (string nme : names){

		string file_name = "../instances/PCSPG-PUCNU/" + nme + ".stp";
		Graph G;
		int _code = STP_reader(file_name, G);
		if (_code != 0){
			cout<<"Error reading file - Code "<<_code<<endl;
			return;
		}

		execute(G, nme);
	}
	cout<<"pcnu  ";
	print_statistics();	
}


void runall(){
	runp();
	runk();
	runc();
	rund();
	rune();
	runh();
	runi640();
}

#endif