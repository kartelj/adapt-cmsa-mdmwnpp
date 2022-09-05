#include "VNS.h"
#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

/// <summary>
/// This is re-implementation of VNS algorithm proposed in A Variable Neighborhood Search Approach for Solving the Multidimensional Multi - Way Number Partitioning Problem
/// by Faria et al.
/// </summary>
/// <param name="problem"></param>
/// <param name="max_seconds"></param>

VNS::VNS(Problem problem, int max_seconds, int r_max)
{
	this->problem = problem;
	this->max_seconds = max_seconds;
	this->r_max = r_max;
	this->start = chrono::steady_clock::now();
}

void VNS::initialize_solution()
{
	s = new int[problem.get_n()];
	sp = new int[problem.get_n()];
	for (int i = 0; i < problem.get_n(); i++)
		s[i] = -1;
	p_card = new int[problem.get_k()];
	p_sum = new double* [problem.get_k()];
	pp_card = new int[problem.get_k()];
	pp_sum = new double* [problem.get_k()];
	for (int i = 0; i < problem.get_k(); i++) {
		p_card[i] = 0;
		p_sum[i] = new double[problem.get_m()];
		pp_sum[i] = new double[problem.get_m()];
		for (int j = 0; j < problem.get_m(); j++) {
			p_sum[i][j] = 0;
			pp_sum[i][j] = 0;
		}
	}
	int l = rand() % problem.get_m();
	for (int i = 0; i < problem.get_n(); i++) {
		int minj = 0;
		int minLj = p_sum[0][l];
		for (int j = 1; j < problem.get_k(); j++) {
			if (p_sum[j][l] < minLj) {
				minLj = p_sum[j][l];
				minj = j;
			}
		}
		s[i] = minj;
		for(int j=0; j<problem.get_m(); j++)
			p_sum[minj][j] += problem.get_v(i, j);
		p_card[minj]++;
	}
	s_fit = problem.fitness(s);
}

void VNS::shake(int r)
{
	// making r independent moves
	for (int q = 0; q < r; q++) {
		int i = -1;
		// finding appropriate i
		while (true) {
			i = rand() % problem.get_n();
			if (pp_card[sp[i]] > 1)
				break; // move allowed
		}
		// now moving to random part
		int p = -1;
		while (true) {
			p = rand() % problem.get_k();
			if (p != sp[i])
				break;
		}
		move(i, p);
	}
}

void VNS::move(int i, int p)
{
	pp_card[sp[i]]--;
	pp_card[p]++;
	for (int j = 0; j < problem.get_m(); j++) {
		pp_sum[sp[i]][j] -= problem.get_v(i, j);
		pp_sum[p][j] += problem.get_v(i, j);
	}
	sp_fit = 0;
	for (int z = 1; z < problem.get_k(); z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem.get_m(); q++) {
				double diff = fabs(pp_sum[z][q] - pp_sum[j][q]); // partial fit. evaluation
				if (diff > sp_fit)
					sp_fit = diff;
			}
	}
	sp[i] = p;
}

void VNS::move_best(int i)
{
	// best is current
	int bestp = sp[i];
	double best_fit = sp_fit;
	for (int p = 0; p < problem.get_k(); p++) {
		if (sp[i] == p)
			continue;
		int old_p = sp[i];
		move(i, p);
		if (sp_fit < best_fit) {
			best_fit = sp_fit;
			bestp = p;
		}
		move(i, old_p);
	}
	if (sp[i] != bestp)
		move(i, bestp);
}

bool VNS::best1()
{
	int besti1 = -1;
	double best = sp_fit;
	for (int i1 = 0; i1 < problem.get_n(); i1++) {
		if (pp_card[sp[i1]] == 1)
			continue; // move not allowed
		int old_p1 = sp[i1];
		move_best(i1);
		if (sp_fit+0.0001 < best) { // to avoid movements because of small acumulation error
			besti1 = i1;
			best = sp_fit;
		}
		move(i1, old_p1); // bringing back to previous state
	}
	if (besti1 != -1) {
		move_best(besti1);
		if (abs(sp_fit - best) > 0.1)
			cout << "Error in best1" << endl;
			//throw new exception("Incorrect calculation in best1.");
		//cout << "best1 improvement" << endl;
		return true;
	}
	return false;
}

bool VNS::best2()
{
	int besti1 = -1;
	int bestp1 = -1;
	int besti2 = -1;
	double best = sp_fit;
	for (int i1 = 0; i1 < problem.get_n(); i1++) {
		if (!time_ok())
			break;
		if (pp_card[sp[i1]] == 1)
			continue; // move not allowed
		int old_p1 = sp[i1];
		for (int p1 = 0; p1 < problem.get_k(); p1++) {
			if (!time_ok())
				break;
			if (sp[i1] == p1)
				continue;
				move(i1,p1);
				for (int i2 = 0; i2 < problem.get_n(); i2++) {
					if (!time_ok())
						break;
					if (i1 == i2 || pp_card[sp[i2]] == 1)
						continue; // move not allowed
					int old_p2 = sp[i2];
					move_best(i2);
					if (sp_fit + 0.0001 < best) {
						bestp1 = p1;
						besti1 = i1;
						besti2 = i2;
						best = sp_fit;
						cout << "New best " << to_string(best) << endl;
					}
					move(i2, old_p2);
				}
				move(i1, old_p1); // bringing back to previous state
		}
	}
	if (besti1 != -1 && besti2 != -1) {
		move(besti1, bestp1);
		move_best(besti2);
		if (abs(sp_fit - best) > 0.1)
			cout << "Error in best2" << endl;
		//throw new exception("Incorrect calculation in best2.");
		cout << "best2 improvement " <<to_string(best) << endl;
		return true;
	}
	return false;
}

bool VNS::best3()
{
	int besti1 = -1;
	int bestp1 = -1;
	int besti2 = -1;
	int bestp2 = -1;
	int besti3 = -1;
	double best = sp_fit;
	for (int i1 = 0; i1 < problem.get_n(); i1++) {
		if (!time_ok())
			break;
		if (pp_card[sp[i1]] == 1)
			continue; // move not allowed
		int old_p1 = sp[i1];
		for (int p1 = 0; p1 < problem.get_k(); p1++) {
			if (!time_ok())
				break;
			if (sp[i1] == p1)
				continue;
			move(i1, p1);
			for (int i2 = 0; i2 < problem.get_n(); i2++) {
				if (!time_ok())
					break;
				if (i1 == i2 || pp_card[sp[i2]] == 1)
					continue; // move not allowed
				int old_p2 = sp[i2];
				for (int p2 = 0; p2 < problem.get_k(); p2++) {
					if (!time_ok())
						break;
					if (sp[i2] == p2)
						continue;
					move(i2, p2);
					for (int i3 = 0; i3 < problem.get_n(); i3++) {
						if (!time_ok())
							break;
						if (i1 == i3 || i2 == i3 || pp_card[sp[i3]] == 1)
							continue; // move not allowed
						int old_p3 = sp[i3];
						move_best(i3);
						if (sp_fit + 0.0001 < best) {
							besti1 = i1;
							bestp1 = p1;
							besti2 = i2;
							bestp2 = p2;
							besti3 = i3;
							best = sp_fit;
						}
						move(i3, old_p3);
					}
					move(i2, old_p2);
				}
			}
			move(i1, old_p1); // bringing back to previous state
		}
	}
	cout << "Finished" << endl;
	if (besti1 != -1 && besti2 != -1 && besti3!=-1) {
		move(besti1, bestp1);
		move(besti2, bestp2);
		move_best(besti3);
		if (abs(sp_fit - best) > 0.1)
			cout << "Error in best3 " << to_string(sp_fit)<<" "<<to_string(best)<< endl;
			//throw new exception("Incorrect calculation in best3.");
		//cout << "best3 improvement" << endl;
		return true;
	}
	return false;
}

void VNS::copy_s_to_sp()
{
	for (int i = 0; i < problem.get_n(); i++)
		sp[i] = s[i];
	for (int i = 0; i < problem.get_k(); i++) {
		pp_card[i] = p_card[i];
		for (int j = 0; j < problem.get_m(); j++)
			pp_sum[i][j] = p_sum[i][j];
	}
	sp_fit = s_fit;
}

void VNS::copy_sp_to_s()
{
	for (int i = 0; i < problem.get_n(); i++)
		s[i] = sp[i];
	for (int i = 0; i < problem.get_k(); i++) {
		p_card[i] = pp_card[i];
		for (int j = 0; j < problem.get_m(); j++)
			p_sum[i][j] = pp_sum[i][j];
	}
	s_fit = sp_fit;
}

bool VNS::time_ok() {
	return chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - this->start).count() < max_seconds;
}

void VNS::run()
{
	int r = 1;
	int it = 0;
	initialize_solution();
	while (time_ok()) {
		copy_s_to_sp();
		shake(r);

		while (time_ok() && best1());

		if (r == 2)
			while(time_ok() && best2());
		else if (r == 3)
			while(time_ok() && best3());
		
		/*
		double sp_fit_check = problem.fitness(sp, k);
		if (abs(sp_fit_check-sp_fit)>0.0001) {
			problem.fitness(sp, k);
			throw new exception("Error in fast fitness evaluation.");
		} */
		if (sp_fit+0.0001 < s_fit) {
			copy_sp_to_s();
			for (int i = 0; i < problem.get_n(); i++)
				s[i] = sp[i];
			r = 1;
			s_fit = sp_fit;
		}
		else {
			r++;
			if (r > r_max)
				r = 1;
		}
		if(it%10==0)
			cout << "Iteration\t" << it <<"\tNeighborhood\t"<<r<< "\tFitness\t" << fixed<< setprecision(4)<< s_fit << endl;
		it++;
	}
	exec_sec = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - this->start).count();
	cout << "Final fitness is " << fixed << setprecision(4) << s_fit << endl;
	cout << "Total time is " << exec_sec << endl;
	cout << "Corresponding partitioning is:" << endl;
	for (int i = 0; i < problem.get_n(); i++)
		cout << s[i] << " ";
	cout << endl;
}

void VNS::output(string path)
{
	ofstream f_out;
	f_out.open(path);
	bool valid = problem.validity_check(s);
	f_out << "value: " << fixed << setprecision(4) << s_fit << endl;
	f_out << "time: " << exec_sec << "\n";
	f_out << "validity: " << valid << endl;
	f_out.close();
}
