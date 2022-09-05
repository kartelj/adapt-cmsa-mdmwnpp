#pragma once
#include <vector>
#include <string>
using namespace std;

class Problem
{
private:
	double **v;
	int n;
	int m;
	int k;

public:
	void load_from_file(string path);
	void load_from_file(string path, int max_n, int max_m, int k);
	double fitness(int sol[]);
	bool validity_check(int* s);
	int dimension();
	int get_m();
	int get_n();
	int get_k();
	double get_v(int i, int j);
};


