#include "Problem.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <set>
using namespace std;

void Problem::load_from_file(string path) {
	string line;
	ifstream f;
	f.open(path);
	if (f.is_open())
	{
		f >> n;
		f >> m;
		f >> k;
		//cout << "Loading n m k from file " << to_string(n) << " " << to_string(m) << " " << to_string(k) << endl;
		v = new double* [n];
		for (int i = 0; i < n; i++) {
			v[i] = new double[m];
			for (int j = 0; j < m; j++)
			{
			        //cout << "Load..."<< i << " " << j <<  endl;
				f >> v[i][j];
			}
		}
		f.close();
	}  
	else cout << "Unable to open file";
 
}

void Problem::load_from_file(string path, int max_n, int max_m, int max_k)
{
	load_from_file(path);
	//cout << "max_n: " << max_n << endl;
	if (max_n > n || max_m > m) {
		cout << "Overriding n or m to values greater than in input files is not allowed." << endl;
		exit(1);
	}
	n = max_n;
	m = max_m;
	k = max_k;
	 /*cout << "Overriding n m k to " << to_string(n) << " " << to_string(m) << " " << to_string(k) << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			cout << v[i][j] << '\t';
		cout << '\n';
	}*/
}

double Problem::fitness(int sol[])
{
	// partition summations
	double** sums = new double* [k];
	for (int i = 0; i < k; i++) {
		sums[i] = new double[m];
		for (int j = 0; j < m; j++)
			sums[i][j] = 0;
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			if (sol[i] != -1)  // sol[i]=-1 has special meaning, so we just skip those (not still set values)
				sums[sol[i]][j] += v[i][j];


	double max_diff = 0;
	// checking all partition summation differences
	for (int i = 1; i < k; i++)
		for (int j = 0; j < i; j++)
			for (int p = 0; p < m; p++) {
				double diff = fabs(sums[i][p] - sums[j][p]);
				if (diff > max_diff)
					max_diff = diff;
			}

	for (int i = 0; i < k; ++i)
		delete[] sums[i];

	return max_diff;
}

bool Problem::validity_check(int* s)
{
	set<int> part_indices;

	for (int i = 0; i < n; i++)
		part_indices.insert(s[i]);

	if (part_indices.size() == k)
		return true;
	else
		return false;
}

int Problem::dimension()
{
	return n;
}

int Problem::get_m()
{
	return m;
}

int Problem::get_n()
{
	return n;
}

int Problem::get_k() {
	return k;
}

double Problem::get_v(int i, int j)
{
	return v[i][j];
}

