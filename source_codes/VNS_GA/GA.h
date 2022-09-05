#pragma once
#include "Problem.h"
#include <tuple>
#include <chrono>
class GA
{
private:
	Problem problem;
	int max_seconds;
	double crossover_rate;
	double mutation_rate;
	int ls3_active;
	long exec_sec;
	int pop_size;
	chrono::steady_clock::time_point start;
	vector<pair<double, int*>> pop;
	void initialize_population();
	int* semi_random_solution();
	double LS(int* sol);
	int select_parent();
	pair<int*, int*> crossover();
	void mutation(int* child);
	double move_fit(int* sp, int i, int p, double** p_sum);
	double swap_fit(int* sp, int i, int j, double** p_sum);
	double swap_fit3(int* sol, int i, int j, int q, double** p_sum);
	bool time_ok();

public:
	GA(Problem problem, int max_seconds, int pop_size, double crossover_rate, double mutation_rate, int ls3_active);
	void run();
	void output(string path);
};

