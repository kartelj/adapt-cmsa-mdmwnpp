#pragma once
#include "Problem.h"
#include <chrono>
class VNS
{
private: 
	Problem problem;
	int r_max;
	int max_seconds;
	long exec_sec;
	chrono::steady_clock::time_point start;
	int *s;
	double s_fit;
	int* p_card; // solution s part cardinality
	double** p_sum; // solution s part sum per coordinate
	int *sp;
	double sp_fit;
	int* pp_card;
	double** pp_sum;
	void initialize_solution();
	void shake(int r);
	void move(int i, int p);
	void move_best(int i);
	bool best1();
	bool best2();
	bool best3();
	void copy_s_to_sp();
	void copy_sp_to_s();
	bool time_ok();

public:
	VNS(Problem problem, int max_seconds, int r_max);
	void run();
	void output(string path);
};

