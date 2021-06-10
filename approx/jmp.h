#ifndef JMP
#define JMP

#include "graph.h"
#include "frac.h"

frac pcst(Graph G);
frac residual_cost(int v, int u);
frac key(int i, int j);
bool insert(int typ, int p, int q);
bool remove(int typ, int p, int q);
void copy(int typ, int p, int q);
frac get_min_key(int typ, int p);
int get_min_ele(int typ, int p);
bool is_in(int typ, int p, int q);
void subcase1b(int p);
void subcase1a(int p, int q);
void init(Graph G);
void iterate();
frac strong_prune(int v, int p, vector<frac> &nw, vector<vector<int>> &adj);
void recover_prune(int v, int p, vector<frac> &nw, vector<vector<int>> &adj);
void prune(Graph G);

#endif