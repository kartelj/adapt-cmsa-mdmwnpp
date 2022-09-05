#include "GA.h"
#include <iostream>
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <unordered_map>

/// <summary>
/// /// This is re-implementation of GA algorithm proposed in A memetic algorithm approach for solving the
/// multidimensional multi - way number partitioning problem by Pop and Matei.
/// </summary>
/// <param name="problem"></param>
/// <param name="max_seconds"></param>
/// <param name="pop_size"></param>
/// <param name="crossover_rate"></param>
/// <param name="mutation_rate"></param>
GA::GA(Problem problem, int max_seconds, int pop_size, double crossover_rate, double mutation_rate, int ls3_active)
{
	this->problem = problem;
	this->max_seconds = max_seconds;
	this->pop_size = pop_size;
	this->crossover_rate = crossover_rate;
	this->mutation_rate = mutation_rate;
	this->ls3_active = ls3_active;
	this->start = chrono::steady_clock::now();
}

void GA::initialize_population()
{
	pop = vector<pair<double, int*>>(pop_size);
	for (int i = 0; i < pop_size; i++) {
		pop[i].second = semi_random_solution();
		pop[i].first = problem.fitness(pop[i].second);
		//cout << to_string(i) << ".\t" << to_string(fit[i]);
		//for (int j = 0; j < problem.get_n(); j++)
		//	cout << pop[i][j] << " ";
		//cout << endl;
	}
}

int* GA::semi_random_solution()
{
	int* sol = new int[problem.get_n()];
	for (int i = 0; i < problem.get_n(); i++)
		sol[i] = -1; // default value
	int q = rand() % (problem.get_n() - 1) + 1;
	// setting genes 1-q randomly
	for (int i = 1; i < q; i++)
		sol[i] = rand() % problem.get_k();
	// others are set to greedily minimize fitness
	for (int i = 0; i < problem.get_n(); i++) {
		if (i >= 1 && i < q)
			continue; // skipping randomly set genes
		int bestj = -1;
		double best = DBL_MAX;
		// checking which part will produce the best fitness
		for (int j = 0; j < problem.get_k(); j++) {
			sol[i] = j;
			double new_fit = problem.fitness(sol);
			if (new_fit < best) {
				best = new_fit;
				bestj = j;
			}
		}
		if (bestj == -1)
			cout << "Error in semi_random_solution()." << endl;
		sol[i] = bestj;
	}
	return sol;
}

double GA::move_fit(int* sol, int i, int p, double** p_sum)
{
	for (int j = 0; j < problem.get_m(); j++) {
		p_sum[sol[i]][j] -= problem.get_v(i, j);
		p_sum[p][j] += problem.get_v(i, j);
	}
	double new_fit = 0;
	for (int z = 1; z < problem.get_k(); z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem.get_m(); q++) {
				double diff = abs(p_sum[z][q] - p_sum[j][q]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;
			}
	}
	for (int j = 0; j < problem.get_m(); j++) {
		p_sum[sol[i]][j] += problem.get_v(i, j);
		p_sum[p][j] -= problem.get_v(i, j);
	}
	return new_fit;
}

double GA::swap_fit(int* sol, int i, int j, double** p_sum)
{
	for (int s = 0; s < problem.get_m(); s++) {
		p_sum[sol[i]][s] -= problem.get_v(i, s);
		p_sum[sol[j]][s] -= problem.get_v(j, s);
		p_sum[sol[i]][s] += problem.get_v(j, s);
		p_sum[sol[j]][s] += problem.get_v(i, s);
	}
	double new_fit = 0;
	for (int z = 1; z < problem.get_k(); z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem.get_m(); q++) {
				double diff = abs(p_sum[z][q] - p_sum[j][q]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;
			}
	}
	for (int s = 0; s < problem.get_m(); s++) {
		p_sum[sol[i]][s] += problem.get_v(i, s);
		p_sum[sol[j]][s] += problem.get_v(j, s);
		p_sum[sol[i]][s] -= problem.get_v(j, s);
		p_sum[sol[j]][s] -= problem.get_v(i, s);
	}
	return new_fit;
}


double GA::swap_fit3(int* sol, int i, int j, int q, double** p_sum)
{
	// i gets what j has 
	// j gets what q has
	// q gets what i has
	for (int s = 0; s < problem.get_m(); s++)
	{
		p_sum[sol[i]][s] -= problem.get_v(i, s);
		p_sum[sol[j]][s] -= problem.get_v(j, s);
		p_sum[sol[q]][s] -= problem.get_v(q, s);
		p_sum[sol[j]][s] += problem.get_v(i, s); // i is moved to partition where is currently j
		p_sum[sol[q]][s] += problem.get_v(j, s); // j is moved to partition where is currently q
		p_sum[sol[i]][s] += problem.get_v(q, s); // q is moved to partition where is currently i
	}
	double new_fit = 0;
	for (int z = 1; z < problem.get_k(); z++) {
		for (int j = 0; j < z; j++)
			for (int r = 0; r < problem.get_m(); r++) {
				double diff = abs(p_sum[z][r] - p_sum[j][r]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;
			}
	}
	//p_sum returns to the same state as before 
	for (int s = 0; s < problem.get_m(); s++) {
		p_sum[sol[i]][s] += problem.get_v(i, s);
		p_sum[sol[j]][s] += problem.get_v(j, s);
		p_sum[sol[q]][s] += problem.get_v(q, s);
		p_sum[sol[j]][s] -= problem.get_v(i, s);
		p_sum[sol[q]][s] -= problem.get_v(j, s);
		p_sum[sol[i]][s] -= problem.get_v(q, s);
	}
	return new_fit;
}

double GA::LS(int* sol)
{
	double fit = problem.fitness(sol);
	double** p_sum = new double* [problem.get_k()];
	for (int i = 0; i < problem.get_k(); i++) {
		p_sum[i] = new double[problem.get_m()];
		for (int j = 0; j < problem.get_m(); j++)
			p_sum[i][j] = 0;
	}
	for (int i = 0; i < problem.get_n(); i++)
		for (int j = 0; j < problem.get_m(); j++)
			p_sum[sol[i]][j] += problem.get_v(i, j);

	// pre-processing counter of elements appearances 
	unordered_map<int, int> maps_count;
	for (int i = 0; i < problem.get_n(); i++)
	{
		int x = sol[i];
		if (maps_count.find(x) == maps_count.end())
			maps_count.insert({ x, 1 });
		else
			maps_count[x]++;
	}

	int impr = 1;
	int n = problem.get_n(); // n = 1
	// It is difficult to understand the explanation in Pop and Matei since 1-change neighbor search says We select randomly an entry from the string representation ... and change its value to some random partition
	// On the other hand complexity is O(n) which means that this is exhaustive LS over all possibilities. 
	// It is probably random in the sense of avoiding position bias (starting from some random vector up to the first improvement). 
	// Similarly for the LS2 and LS3...
	while (impr && time_ok()) {
		impr = 0;
		// LS1
		int i_r = rand() % problem.get_n();
		for (int ix = 0; ix < n; ix++) {
			int i = (i_r + ix) % problem.get_n(); // to avoid positional bias
			// skipping if sol[i] partition is occuring once, in order to avoid going to infeasible solution
			if (maps_count[sol[i]] == 1)
				continue;
			int p = sol[i];
			do {
				p = rand() % problem.get_k();
			} while (p == sol[i]);
			int oldp = sol[i];
			double new_fit = move_fit(sol, i, p, p_sum);
			if (new_fit < fit) {
				for (int j = 0; j < problem.get_m(); j++) {
					p_sum[sol[i]][j] -= problem.get_v(i, j);
					p_sum[p][j] += problem.get_v(i, j);
				}
				// adjusting counts
				maps_count[sol[i]]--;
				maps_count[p]++;
				sol[i] = p;
				fit = new_fit;
				impr = 1;
				//double control_fit = problem.fitness(sol, k);
				//if (fabs(control_fit - fit) > 0.00001)
				//	throw new exception("Incorrect partial fitness function.");
				//cout << "LS1 " << to_string(new_fit) << endl;
				break;
			}
		}
	}
	impr = 1;
	while (impr && time_ok()) {
		impr = 0;
		// LS2
		int ir = rand() % problem.get_n();
		for (int ix = 0; ix < n; ix++) {
			int i = (ir + ix) % problem.get_n(); // to avoid positional bias
			int jr = rand() % problem.get_n();
			for (int jx = 0; jx < n; jx++) {
				int j = (jr + jx) % problem.get_n();
				if (i >= j || sol[i] == sol[j])
					continue;
				double new_fit = swap_fit(sol, i, j, p_sum);
				if (fit - new_fit > 0.1) {
					for (int s = 0; s < problem.get_m(); s++) {
						p_sum[sol[i]][s] -= problem.get_v(i, s);
						p_sum[sol[j]][s] -= problem.get_v(j, s);
						p_sum[sol[i]][s] += problem.get_v(j, s);
						p_sum[sol[j]][s] += problem.get_v(i, s);
					}
					int pi = sol[i];
					sol[i] = sol[j];
					sol[j] = pi;
					fit = new_fit;
					impr = 1;
					//	double control_fit = problem.fitness(sol, k);
					//	if (fabs(control_fit - fit) > 0.001)
					//		throw new exception("Incorrect partial fitness function.");
						//cout << "LS2 " << to_string(new_fit) << endl;
					break;
				}
			}
			if (impr)
				break;
		}
	}
	//LS3
	
	if (ls3_active) {
		impr = 1;
		while (impr && time_ok()) {
			impr = 0;

			int ir = rand() % problem.get_n();

			for (int ix = 0; ix < n; ix++) {

				int i = (ir + ix) % problem.get_n(); // to avoid positional bias
				int jr = rand() % problem.get_n();

				for (int jx = 0; jx < n; jx++) {
					int j = (jr + jx) % problem.get_n();
					if (i >= j || sol[i] == sol[j])
						continue;
					int jr = rand() % problem.get_n();
					for (int qx = 0; qx < n; qx++) { // resove a bias, start from a random position jr
							//da li treba randomiziran  q
						int q = (jr + qx) % problem.get_n();

						if (sol[q] == sol[i] || sol[q] == sol[j])
							continue;

						if (i == 10 && j == 12 && q == 45 && ir == 21 && jr==5)
							cout << "bla" << endl;
						double new_fit = swap_fit3(sol, i, j, q, p_sum);
						if (fit - new_fit > 0.1) {
							// i gets what j has 
							// j gets what q has
							// q gets what i has
							for (int s = 0; s < problem.get_m(); s++)
							{
								p_sum[sol[i]][s] -= problem.get_v(i, s);
								p_sum[sol[j]][s] -= problem.get_v(j, s);
								p_sum[sol[q]][s] -= problem.get_v(q, s);
								p_sum[sol[j]][s] += problem.get_v(i, s); // i is moved to partition where is currently j
								p_sum[sol[q]][s] += problem.get_v(j, s); // j is moved to partition where is currently q
								p_sum[sol[i]][s] += problem.get_v(q, s); // q is moved to partition where is currently i
							}
							//best_fit = new_fit;
							int pi = sol[i];
							int pj = sol[j];
							int pq = sol[q];
							sol[i] = pj;
							sol[j] = pq;
							sol[q] = pi;
							fit = new_fit;
							impr = 1;
							double control_fit = problem.fitness(sol);
							if (fabs(control_fit - fit) > 0.001)
								throw out_of_range("Incorrect partial fitness function.");
							break;
						}
					}
					if (impr)
						break;
				}
			}
		}
	}

	for (int i = 0; i < problem.get_k(); i++)
		delete[] p_sum[i];

	return fit;
}

int GA::select_parent()
{
	int i = rand() % pop_size;
	int j = rand() % pop_size;
	if (pop[i].first < pop[j].first)
		return i;
	else
		return j;
}

pair<int*, int*> GA::crossover()
{
	int p1i = select_parent();
	int p2i = select_parent();
	int* p1 = pop[p1i].second;
	int* p2 = pop[p2i].second;
	double crossover_rv = rand() * 1.0 / RAND_MAX;
	// with chance 1-crossover_rate just return the parents
	if (crossover_rv > crossover_rate)
		return make_pair(p1, p2);
	// otherwise performe 1-point crossover
	int* c1 = new int[problem.get_n()];
	int* c2 = new int[problem.get_n()];
	int pt = rand() % problem.get_n();
	for (int i = 0; i < pt; i++) {
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
	for (int i = pt; i < problem.get_n(); i++) {
		c1[i] = p2[i];
		c2[i] = p1[i];
	}
	return make_pair(c1, c2);
}

void GA::mutation(int* child)
{
	for (int i = 0; i < problem.get_n(); i++)
		if (rand() * 1.0 / RAND_MAX < 0.1) {
			int newp = child[i];
			do {
				newp = rand() % problem.get_k();
			} while (newp == child[i]);
			child[i] = newp;
		}
}

bool GA::time_ok() {
	return chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - this->start).count() < max_seconds;
}

void GA::run()
{
	initialize_population();
	for (int i = 0; i < pop_size; i++)
		pop[i].first = LS(pop[i].second);

	int interm_pop_size = pop_size;
	int g = 0;
	while (time_ok()) {
		// now oversample 10x active population to form intermediate population
		vector<pair<double, int*>> intermediate_pop = vector<pair<double, int*>>(interm_pop_size);
		// elitism -- just copy the best two from the current population
		sort(pop.begin(), pop.end());
		intermediate_pop[0].first = pop[0].first;
		intermediate_pop[0].second = pop[0].second;
		intermediate_pop[1].first = pop[1].first;
		intermediate_pop[1].second = pop[1].second;

		for (int i = 1; i < interm_pop_size / 2; i++) {
			pair<int*, int*> children = crossover();
			mutation(children.first);
			mutation(children.second);
			double fit1 = LS(children.first);
			intermediate_pop[2 * i] = make_pair(fit1, children.first);
			double fit2 = LS(children.second);
			intermediate_pop[2 * i + 1] = make_pair(fit2, children.second);
		}
		// sorting w.r.t. fitness function
		sort(intermediate_pop.begin(), intermediate_pop.end());
		// taking only the best pop_size for the next generation
		for (int i = 0; i < pop_size; i++) {
			//delete pop[i].second;
			pop[i].second = intermediate_pop[i].second;
			pop[i].first = intermediate_pop[i].first;
		}
		/*for (int i = 0; i < 1; i++)
			cout << to_string(g) << ". " << to_string(pop[i].first) << endl;*/
		//cout << "----------------------------------" << endl;
		g++;
		// free up memory
		for (int i = pop_size; i < intermediate_pop.size(); ++i)
			delete[] intermediate_pop[i].second;

	}
	exec_sec = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - this->start).count();
}

void GA::output(string path)
{
	ofstream f_out;
	f_out.open(path);
	//if(pop.size() > 0)
	//{
	   bool valid = problem.validity_check(pop[0].second);
	   cout << "value: " << fixed << setprecision(4) << pop[0].first << endl;
	 f_out << "value: "  << fixed << setprecision(4) << pop[0].first << endl;
	 f_out << "time: " << exec_sec << "\n";
	 f_out << "validity: " << valid << endl;
	 f_out.close();
	//}else{
	   
	//   cout <<"value: "<< 10000000000000 << endl;
	
	//}##### sesija 0: GA za n=100
	// ##### sesija 1: GA za n=50

}
