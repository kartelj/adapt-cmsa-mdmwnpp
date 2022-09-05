#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Timer.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include<windows.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#ifdef _WIN32
#include "ilcplex/ilocplex.h"
#include <string>
#include <iostream>
#include <fstream>
#include <float.h>
#else
//#include <ilcp/cp.h>
// the following "include" is necessary for the correct working/compilation of CPLEX. You should adapt this path to your installation of CPLEX
// #include "/home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/cplex/include/ilcplex/ilocplex.h" 
//"/home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/include/ilcplex/ilocplex.h"
#include "/home/marko/Desktop/CPLEX_Studio127/cplex/include/ilcplex/ilocplex.h"
// "/home1/share/ILOG/cplex-12.7.0/cplex/include/ilcplex/ilocplex.h"
#endif

#define wmat(i,j) problem->w[(i)*problem->m + j]

#define INF 1000000000000000.0 
#define MAXS 300



// using CP; 
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;

/** data structure for Problem **/
typedef struct {
	int n, m, k, d, alg;		/* broj elemenata i dimenzija */
	double* w, * sum;         	/* matrix weights i array of sums of weights */
	double* y;             		/* solutino */
	double fv;              	/* value of obj function */
} problem_struct;



int usl, prom, n;
int nprom, nusl, zero, nprommax, nuslmax;
FILE* ul, * izl, * mst;
//CPXFILEptr log ;
char* tip, ** imepr;
int* ind;
//std::string imeul; 
char imeul[MAXS], imeizl[MAXS], imelog[MAXS], imemst[MAXS], imepr2[MAXS];
clock_t rs, re;
double vreme, rv;

double cur_time = 0.0;
int iter_move = 1;
double best_rand_sol = INF;
/* CMSA parameters */
int cmsa_cplex_time = 10;
int n_a = 50;
int age_max = 5;
int cmsa_milp = 0;  //COAM by default
int cmsa_greedy = 2; //KMEANS by default

vector<int> s_bsf; // store  best solution
double obj_best = INF; //best obj in CMSA
int indeks = 0; // for multiple runs of an instance 
int seed = 1;
/** end of CMSA parameters **/

/** parameters of Adapted-CMSA **/

double alpha_LB = 0.35;
double alpha_UB = 0.97;
double alpha_red = 0.05;
double t_prop = 0.13;
double ls_number = 2; // how many different LS is used
double prefer_s_bsf = 1;

/** end Adapt-CMSA parameters definition **/

// VNS CMSA parameter
double move_prob = 0.5;
int min_p = 2;
int max_p = 50;

Timer timer;
problem_struct* problem;


/** Enum types for algorithms **/
enum MILPS {

	COAM = 0,
	FARIA = 1

};

enum HEURISTICS
{

	CMSAH = 4,
	ADAPT_CMSAH = 5,
	DEEP_CMSAH = 6
};

enum GREEDY
{
	KMEANG = 2,
	MDRGHG = 3
};
/** end of Enum types **/

string outPath = "/home/marko/Desktop/Papers/MTWNPP/MTWNPP/CMSA/code/";

ILOSOLVECALLBACK2(abortCallback, IloCplex::Aborter&, abo, double&, curbest) {

	if (hasIncumbent()) {

		IloNum nv = getIncumbentObjValue();

		if (curbest > double(nv))

			abo.abort(); // this assumes you are minimizing the objective function
	}
}


/* input data  */

void ulazpod(void)
{
	int i, j, k, l;
	double pom;

	ul = fopen(imeul, "rt");
	if (ul == NULL)
	{
		printf("Ulazna datoteka nije otvorena!\n");
		exit(0);
	}

	if (problem->n != NULL and problem->m != NULL and problem->k != NULL)
	{
		fscanf(ul, "%d %d", &k, &l);
		if (k < problem->n || l < problem->m)
		{
			printf("U datoteci nema dovoljno brojeva!\n");
			exit(0);
		}
	}
	else {

		fscanf(ul, "%d %d %d", &problem->n, &problem->m, &problem->k);
		cout << "n: " << problem->n << " " << problem->m << " " << problem->k << endl;
		k = problem->n;
		l = problem->m;
	}

	problem->w = (double*)malloc(problem->n * problem->m * sizeof(double));
	problem->y = (double*)malloc(problem->n * sizeof(double));
	problem->sum = (double*)malloc(problem->m * sizeof(double));

	if (problem->w == NULL || problem->y == NULL || problem->sum == NULL)
	{
		printf("Dinamicke promenljive nisu dobro alocirane!\n");
		exit(0);
	}

	for (i = 0; i < problem->n; i++)
	{
		for (j = 0; j < problem->m; j++)
			fscanf(ul, "%lf", &wmat(i, j));
		for (; j < l; j++)
			fscanf(ul, "%lf", &pom);
	}

	for (j = 0; j < problem->m; j++)
	{
		problem->sum[j] = 0;
		for (i = 0; i < problem->n; i++)
			problem->sum[j] += wmat(i, j);
	}
	/*for (i = 0; i < problem->n; i++)
	{
		for (j = 0; j < problem->m; j++)
			cout << wmat(i, j) << "\t ";
		cout << endl;
	}*/
	fclose(ul);
}
/* ulazpod */






double max_minus_min(vector<vector<double>>& sum, int l)
{
	//cout << "max min " << endl;
	double max_v = -10000000.0; double min_v = 1000000000.0;

	for (int j = 0; j < problem->k; ++j)
	{

		if (max_v < sum[j][l])
			max_v = sum[j][l];
	}
	//cout << "max_v : " << max_v << endl;

	for (int j = 0; j < problem->k; ++j)
	{

		if (min_v > sum[j][l])
			min_v = sum[j][l];
	}
	return max_v - min_v;

}


double objective(vector<int>& solution, int k = -1)
{

	if (solution.empty())
		return INF;

	double value = 0.0;
	vector<vector<double>> sum;  // sum[i][dx]: sum of coordinate dx from {0,...,m-1} of vectors in partition i 

	int kval = k;

	if (kval == -1) // if alg == 3 (KMeans based) ==> problem->k is replaced to k (allowed number of partitions)
		kval = problem->k;

	for (int i = 0; i < kval; ++i)
	{
		vector<double> sum_i;
		for (int j = 0; j < problem->m; ++j)
			sum_i.push_back(0.0);

		sum.push_back(sum_i);
	}

	int index = 0;
	for (int px : solution)
	{

		for (int j = 0; j < problem->m; ++j)
			sum[px][j] += wmat(index, j);

		index++; // next vector
	}
	// calculate min max value:
	//cout << "min max calc" << endl;
	double obj = -INF;

	for (int j = 0; j < problem->m; ++j) // through coordinates:
	{
		double maxmin = max_minus_min(sum, j);
		if (obj < maxmin)
			obj = maxmin;
	}

	return obj;

}


bool validityCheck(vector<int>& s)
{
	set<int> part_indices;

	for (auto& index : s)
		part_indices.insert(index);

	if (part_indices.size() == problem->k)
		return true;
	else
		return false;

}

/* Kojic, m=2 CPLEX  */
void cplex_init_kojic()
{
	/*for(int i = 0; i < problem->m; ++i)
		cout << "sum[ " << i << " ] " << problem->sum[ i ] << "\n" ;
	*/
	IloEnv env;
	IloModel model(env);

	IloNumVar z(env);
	IloNumVarArray y(env);
	//cout << "Vars initialized " << endl; 
	for (int i = 0; i < problem->n; ++i)
		y.add(IloNumVar(env, 0, 1, ILOINT));

	//cout << "Constriants " << endl; 
	/* uslov (1) */
	for (int j = 0; j < problem->m; j++) // iterate through dimensions: 
	{
		IloExpr e(env);
		for (int ix = 0; ix < problem->n; ++ix)
			e += wmat(ix, j) * y[ix];

		e -= 0.5 * z;
		model.add(e <= 0.5 * problem->sum[j]);
	}
	/** uslov (2) **/
	for (int j = 0; j < problem->m; j++) // iterate through dimensions: 
	{
		IloExpr e(env);
		for (int ix = 0; ix < problem->n; ++ix)
			e += wmat(ix, j) * y[ix];
		e += 0.5 * z;
		model.add(e >= 0.5 * problem->sum[j]);
	}
	/**  uslov 3 iz Faria (2017) **/
	/*IloExpr expr(env);
	for(int j = 0; j < problem->n; j++) // iterate through dimensions:
	{
		expr += y[ j ];
	}*/
	//model.add(expr >= 1); 
	model.add(IloMinimize(env, z)); // max function
	IloCplex cplex(model);
	// cplex configuration:
	int time_limit = vreme;
	// pass the time limit to CPLEX
	cplex.setParam(IloCplex::TiLim, time_limit);

	// the following two parameters should always be set in the way as shown
	cplex.setParam(IloCplex::NodeFileInd, 2);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::EpGap, 0.0000000001);
	cplex.setParam(IloCplex::EpAGap, 0.00001);

#ifdef _WIN32
	IloNum lastObjVal = DBL_MAX;
#else
	IloNum lastObjVal = std::numeric_limits<double>::max();
#endif
	//cplex.exportModel("lpex1.lp");
	// tell CPLEX to make use of the function 'loggingCallback' for writing out information to the screen
	//cplex.use(loggingCallback(env, timer,/* times, results, gaps, iter,*/ lastObjVal));
	// rjesavanje modela:
	cplex.solve();

	if (cplex.getStatus() == IloAlgorithm::Optimal or cplex.getStatus() == IloAlgorithm::Feasible)
	{
		/*if(cplex.getStatus() == IloAlgorithm::Optimal)
		   cout << "CPLEX finds optimal" << endl;
		else
		   cout << "CPLEX finds feasible solution" << endl;*/

		double lastVal = double(cplex.getObjValue());
		// print the objective point
		cout << "nodes/vertices in the solution: {" << endl;
		bool first = true;
		for (int i = 0; i < problem->n; ++i) {
			IloNum xval = cplex.getValue(y[i]);

			if (xval > 0.9) {
				cout << "vector " << i << " included into partitioning 1 " << endl;
			}
		}
		cout << "}" << endl;
		cout << "value: " << lastVal << endl;
		cout << "gap: " << double(cplex.getMIPRelativeGap()) << endl;
		double end_time = timer.elapsed_time(Timer::VIRTUAL);
		cout << "time: " << (end_time - cur_time) << "\n";

		ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m) +
			"_" + std::to_string(problem->k) + "_" + std::to_string(problem->alg) + ".out");
		// write  into a file
		myfileOut << "value: " << lastVal << endl;
		myfileOut << "time: " << (end_time - cur_time) << "\n";
		myfileOut << "gap: " << double(cplex.getMIPRelativeGap()) << endl;
	}

	env.end();
}


// Apdapt Faria model bin. vars ==> x_{ji} ==> x_{ij}
vector<int> cplex_init_faria(set<pair<int, int>>& C_prime, double UB = INF)
{

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	IloCplex::Aborter abo(env);
	// set up aborter
	if (problem->alg == CMSAH) // or problem->alg == 5)
	{
		cplex.use(abo);  // here you tell the IloCplex object to use this aborter
		cplex.use(abortCallback(env, abo, obj_best));  // here you tell the IloCplex object to use "your" aborter function from above.  
	}
	IloNumVar y(env);

	vector<int> solution;

	NumVarMatrix x(env, problem->n);
	for (int i = 0; i < problem->n; i++) {
		x[i] = IloNumVarArray(env, problem->k);
		for (int k = 0; k < problem->k; ++k)
			x[i][k] = IloNumVar(env, 0.0, 1.0, ILOINT);
	}
	//cout << "Vars initialized in Faria model... " << endl; 
	// cout << "Now constriants... " << endl; 

	 // if CMSA has been executed ==> form a subinstance to solve 
	if (!C_prime.empty())
	{

		for (int i = 0; i < problem->n; i++) {
			for (int j = 0; j < problem->k; j++) {

				pair<int, int> c = make_pair(i, j);
				if (C_prime.find(c) == C_prime.end())
				{
					model.add(x[i][j] == 0);
				}
			}
		}
	}

	/** 11c **/
	for (int j = 0; j < problem->k; ++j)
	{
		IloExpr e(env);
		for (int i = 0; i < problem->n; ++i)
		{
			e += x[i][j];
		}
		model.add(e >= 1);
	}
	/** 11b  **/
	for (int i = 0; i < problem->n; ++i)
	{
		IloExpr e(env);
		for (int j = 0; j < problem->k; ++j)
		{
			e += x[i][j];
		}
		model.add(e == 1);

	}
	/** 11d **/
	for (int j1 = 0; j1 < problem->k; ++j1)
	{
		for (int j2 = j1 + 1; j2 < problem->k; ++j2)
		{
			for (int l = 0; l < problem->m; ++l)
			{
				IloExpr e_ijl(env);
				for (int i = 0; i < problem->n; ++i)
				{
					e_ijl += wmat(i, l) * (x[i][j1] - x[i][j2]);
				}
				model.add(e_ijl <= y);
			}
		}
	}
	/**  11e: symetric to 11d **/
	for (int j1 = 0; j1 < problem->k; ++j1)
	{
		for (int j2 = j1 + 1; j2 < problem->k; ++j2)
		{
			for (int l = 0; l < problem->m; ++l)
			{
				IloExpr e_ijl(env);
				for (int i = 0; i < problem->n; ++i)
				{
					e_ijl += wmat(i, l) * (x[i][j2] - x[i][j1]);
				}
				model.add(e_ijl <= y);
			}
		}
	}

	/** add the objective **/
	model.add(IloMinimize(env, y)); // 11a function

	/** cut-off **/
	if (UB != INF)
		model.add(y <= UB);  // apply a cut-off based on UB

	  /** setUp Cplex parameters and solve the  model **/
	int time_limit = vreme;
	if (problem->alg == CMSAH or problem->alg == ADAPT_CMSAH or problem->alg == DEEP_CMSAH) // a subinstance has been solved (CMSA is executed):
	{

		double time_stop = timer.elapsed_time(Timer::VIRTUAL);
		if (vreme - time_stop < cmsa_cplex_time)
			time_limit = vreme - time_stop;
		else
			time_limit = cmsa_cplex_time;
	}
	if (time_limit < 0.9)
	{
		time_limit = 1.0; // time exceeded 
		// cmsa_cplex_time = 0;  
	}
	// pass the time limit to CPLEX
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, time_limit);
	// the following two parameters should always be set in the way as shown

	cplex.setParam(IloCplex::NodeFileInd, 2);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::EpGap, 0.00001);
	cplex.setParam(IloCplex::EpAGap, 0.00001);

#ifdef _WIN32
	IloNum lastObjVal = DBL_MAX;
#else
	IloNum lastObjVal = std::numeric_limits<double>::max();
#endif
	cplex.solve();
	/** check if a model returs useful solutions **/
	if (cplex.getStatus() == IloAlgorithm::Optimal or cplex.getStatus() == IloAlgorithm::Feasible)
	{
		/*if(cplex.getStatus() == IloAlgorithm::Optimal)
		   cout << "CPLEX finds optimal" << endl;
		else
		   cout << "CPLEX finds feasible solution" << endl;*/

		double lastVal = double(cplex.getObjValue());
		// print the objective point

		cout << "nodes/vertices in the solution: {" << endl;
		bool first = true;
		for (int i = 0; i < problem->n; ++i) {
			for (int j = 0; j < problem->k; ++j) {

				IloNum xval = cplex.getValue(x[i][j]);
				if (xval > 0.9) {
					//cout << "Vector " << i << " included into partitioning "  << j << endl; 
					solution.push_back(j);
				}
			}
		}
		cout << "}" << endl;
		cout << "value: " << lastVal << endl;
		cout << "gap: " << double(cplex.getMIPRelativeGap()) << endl;
		double end_time = timer.elapsed_time(Timer::VIRTUAL);
		cout << "time: " << (end_time - cur_time) << "\n";
		if (problem->alg == FARIA) // only if the pure MILP is run, write the solution in an out file
		{
			ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m) +
				"_" + std::to_string(problem->k) + "_" + std::to_string(problem->alg) + ".out");
			// write  into a file
			myfileOut << "value: " << std::fixed << objective(solution) << endl; // resolve discrepancy in the final objective value
			myfileOut << "time: " << std::fixed << (end_time - cur_time) << "\n";
			myfileOut << "gap: " << std::fixed << double(cplex.getMIPRelativeGap()) << endl;
		}
	}
	env.end();


	return solution;
}



// find argmin values of L
int argmin(vector<vector<double>>& L, vector<int>& D)
{
	int j_star = -1;
	double val = INF;
	for (int j = 0; j < problem->k; ++j) // j in {0,...,k-1}
	{
		for (int dx = 0; dx < problem->d; ++dx)
		{
			if (val > L[j][dx])
			{
				val = L[j][dx];
				j_star = j; //argmin
			}
		}
	}

	return j_star;
}


/** COAM model, CPLEX **/
vector<int>  cplex_COAM(set<pair<int, int>>& C_prime, double UB = INF)
{

	IloEnv env;
	IloModel model(env);

	IloCplex::Aborter abo(env);  // this initializes an aborter

	IloCplex cplex(model);

	if (problem->alg == CMSAH) //or problem->alg == 5)
	{
		cplex.use(abo);  // here you tell the IloCplex object to use this aborter
		cplex.use(abortCallback(env, abo, obj_best));  // here you tell the IloCplex object to use "your" aborter function from above.  
	}

	vector<int> solution; // save solution in the form of a vector
	//cout << "Vars initialisation in COAM model " << endl; 
	// x_ij \in {0, 1}
	NumVarMatrix x(env, problem->n);
	for (int i = 0; i < problem->n; i++) {
		x[i] = IloNumVarArray(env, problem->k);
		for (int j = 0; j < problem->k; ++j)
			x[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT);
	}

	// if CMSA has been executed ==> form a subinstance to solve 
	if (!C_prime.empty())
	{
		int added = 0;
		for (int i = 0; i < problem->n; i++) {
			for (int j = 0; j < problem->k; j++) {

				pair<int, int> c = make_pair(i, j);
				if (C_prime.find(c) == C_prime.end())
				{
					model.add(x[i][j] == 0); added++;
				}
			}
		}
		//cout << "added ..." << added << " constraint " << UB << endl;
	}

	// variables y_l
	IloNumVarArray y(env, problem->m);
	for (int l = 0; l < problem->m; ++l)
		y[l] = IloNumVar(env, -IloInfinity, IloInfinity); // y_l \in R

	// variables z_l
	IloNumVarArray z(env, problem->m);
	for (int l = 0; l < problem->m; ++l)
		z[l] = IloNumVar(env, -IloInfinity, IloInfinity); // z_l \in R

	// minimize r:
	IloNumVar r(env, 0, IloInfinity);
	model.add(IloMinimize(env, r)); // (1)

	if (UB != INF) //cut-off the search space by appropriate UB
		model.add(r <= UB);
	//cout << "Constriants " << endl;  
	// Constr. (2) 
	for (int i = 0; i < problem->n; ++i)
	{
		IloExpr e(env);
		for (int j = 0; j < problem->k; ++j)
		{
			e += x[i][j];
		}
		model.add(e == 1);
	}
	// Constr. (3) 
	for (int j = 0; j < problem->k; ++j)
	{
		for (int l = 0; l < problem->m; ++l)
		{
			IloExpr e(env);

			for (int i = 0; i < problem->n; ++i)
				e += wmat(i, l) * x[i][j];

			model.add(e <= y[l]);
		}
	}
	// Constr. (4)
	for (int j = 0; j < problem->k; ++j)
	{
		for (int l = 0; l < problem->m; ++l)
		{
			IloExpr e(env);

			for (int i = 0; i < problem->n; ++i)
				e += wmat(i, l) * x[i][j];

			model.add(e >= z[l]);
		}
	}
	// Constr. (5)
	for (int l = 0; l < problem->m; ++l)
	{

		model.add(y[l] - z[l] <= r);

	}
	//Faria et al. model --> (included by desire) 
	for (int j = 0; j < problem->k; ++j)
	{
		IloExpr e(env);
		for (int i = 0; i < problem->n; ++i)
			e += x[i][j];
		model.add(e >= 1);
	}


	// Cplex settings 
	int time_limit = vreme;
	if (problem->alg == CMSAH or problem->alg == ADAPT_CMSAH or problem->alg == DEEP_CMSAH  )  // a subinstance has been solved (CMSA is executed):
	{

		double time_stop = timer.elapsed_time(Timer::VIRTUAL);
		if (vreme - time_stop < cmsa_cplex_time)
			time_limit = vreme - time_stop;
		else
			time_limit = cmsa_cplex_time;
	}
	if (time_limit < 0.9)
	{
		time_limit = 1; // time exceeded (tolerance of 0.1 s)
	   //cmsa_cplex_time = 0.0; 
	}
	// pass the time limit to CPLEX
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, time_limit);//CPX_PARAM_EPINT
	// the following two parameters should always be set in the way as shown
	cplex.setParam(IloCplex::NodeFileInd, 2);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::EpGap, 0.000000000000001);//EpAGap
	cplex.setParam(IloCplex::EpAGap, 0.000000000000001);
	//cplex.setParam(IloCplex::EpInt, 0.000000000000001);//EpAGap
	//cplex.setParam(IloCplex::NumericalEmphasis, 1);//EpAGap

#ifdef _WIN32
	IloNum lastObjVal = DBL_MAX;
#else
	IloNum lastObjVal = std::numeric_limits<double>::max();
#endif
	cplex.solve();
	/** check if a model returs useful solutions **/
	if (cplex.getStatus() == IloAlgorithm::Optimal or cplex.getStatus() == IloAlgorithm::Feasible)
	{
		/*if(cplex.getStatus() == IloAlgorithm::Optimal)
		   cout << "CPLEX finds optimal" << endl;
		else
		   cout << "CPLEX finds feasible solution" << endl;*/

		double lastVal = double(cplex.getObjValue());
		// print the objective point
		//cout << "nodes/vertices in the solution: {" <<endl;
		bool first = true;
		for (int i = 0; i < problem->n; ++i) {
			for (int j = 0; j < problem->k; ++j) {

				IloNum xval = cplex.getValue(x[i][j]);

				if (xval > 0.9) {
					//cout << "Vector " << i << " included into partitioning "   << j << endl; 
					solution.push_back(j);
				}
			}
		}
		//cout << "}" << endl;
		//cout << "value: " << lastVal << endl;
		//cout << "gap: " << double(cplex.getMIPRelativeGap()) << endl;
		double end_time = timer.elapsed_time(Timer::VIRTUAL);
		//cout << "time: " << (end_time - cur_time ) <<"\n";
		//assert(lastVal == objective(solution) );
		//cout << "Validation of the solution: " << objective(solution) << endl; 

		if (problem->alg == COAM) // only if the pure MILP is run, write the solution in an out file
		{
			ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m) +
				"_" + std::to_string(problem->k) + "_" + std::to_string(problem->alg) + ".out");
			// write  into a file
			myfileOut << "value: " << std::fixed << objective(solution) << endl;
			myfileOut << "time: " << std::fixed << (end_time - cur_time) << "\n";
			myfileOut << "gap: " << std::fixed << double(cplex.getMIPRelativeGap()) << endl;
		}

	}
	env.end();

	return solution;
}


/** Greedy heuristic: MDRGH
	d: input parameter ==> number of dimensions used to guide search
**/

// rearrange the solution according to a partial ordering relation (the smallest index of vector in partition i is smaller than that of partition j, for i < j).   
void rearrange(vector<int>& solution)
{

	vector<int> rearrange(problem->n, -1);
	int index = 0; // indexing partitions from 0 to k
	for (int i = 0; i < problem->n; ++i)
	{
		int px = solution[i];
		if (rearrange[i] == -1)
		{
			rearrange[i] = index;
			for (int j = i + 1; j < problem->n; ++j)
				if (solution[j] == px)
					rearrange[j] = index;

			index++;
		}
	}
	index = 0;
	for (auto x : rearrange)
	{
		solution[index] = x;
		index++;
	}
}


vector<int>  MDRGH(int d)
{
	cout << "MDRGH executed..." << endl;

	vector<vector<double>> L;
	vector<int> solution(problem->n, -1);

	for (int j = 0; j < problem->k; ++j)
	{
		vector<double> L_j;
		for (int dx = 0; dx < problem->d; ++dx)
			L_j.push_back(0.0);

		L.push_back(L_j);
	}

	vector<int> D; // store random dymensions utilized to guide search 
	//srand(time(NULL));

	for (int dx = 0; dx < d; ++dx)
	{
		int dim = (rand() % problem->m); // rand 

		if (std::find(D.begin(), D.end(), dim) == D.end())
			D.push_back(dim);
		else
			dx--;
	}
	// iterate
	vector<int> shuffle_arr;
	for (int i = 0; i < problem->n; ++i)
		shuffle_arr.push_back(i);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	shuffle(shuffle_arr.begin(), shuffle_arr.end(), std::default_random_engine(seed));

	for (int i = 0; i < problem->n; ++i) // shuffle 
	{
		//argmin L-values  
		int j_star = argmin(L, D);
		solution[shuffle_arr[i]] = j_star;
		// update L-values

		for (int dx = 0; dx < d; ++dx) // problem
			L[j_star][dx] = L[j_star][dx] + wmat(shuffle_arr[i], D[dx]); // i ==>  shuffle [ i ] (shuffled) 

	}

	double obj = objective(solution);

	//cout << "Rearranged solution: " << endl;
	rearrange(solution);

	if (problem->alg == MDRGHG)
	{
		double end_time = timer.elapsed_time(Timer::VIRTUAL);
		cout << "time: " << (end_time - cur_time) << "\n";
		cout << " obj: " << objective(solution) << endl;
		// store the solution: 
		ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m) +
			"_" + std::to_string(problem->k) + "_" + std::to_string(problem->d) + "_" + std::to_string(problem->alg) + ".out");
		// write  into a file

		myfileOut << "obj: " << objective(solution) << endl;
		myfileOut << "value: " << std::fixed << obj << endl;
		myfileOut << "time: " << (end_time - cur_time) << "\n";
		myfileOut.close();
	}
	//cout << "obj: " << objective(solution) << endl;
	return solution;
}

vector<int> RandomizedGenerator()
{
	//cout << "Randomized procedure... " << endl;
	// start with n-partitioning --> each vector in its own partition
	vector<int> solution(problem->n, -1);

	int kx = 0;
	while (kx < problem->k) // ensure feasibility 
	{
		int k = rand() % problem->n;
		while (solution[k] != -1)
		{
			k = rand() % problem->n; cout << "k: " << k << endl;
		}
		solution[k] = kx;

		kx++;
	}
	for (auto i = 0; i < problem->n; ++i)  // geenrate random solution
	{
		if (solution[i] == -1)
		{

			int part_r = rand() % problem->k;
			solution[i] = part_r;
		}
	}
	//double obj = LSbest(solution);
	rearrange(solution);
	return solution;

}

vector<int> KMeansHeuristic(int num_move = 1)
{
	cout << "KMeans based heuristic " << endl;
	// start with n-partitioning --> each vector in its own partition
	vector<int> solution;

	for (auto i = 0; i < problem->n; ++i)
		solution.push_back(i);

	//merge partitioning ==> two vectors into one randomly 
	int count_merge = 0;

	while (true)
	{

		// Neighborhood ==> two random partitioning and merge them into one (try multiple timesand choose the best acc to obj value) before move to merging)...
		int move = 0;
		double best_obj_move = INF;
		vector<int> solution_prime;
		int move_v1 = -1;
		int move_v2 = -1;

		while (move < num_move)
		{

			for (auto& x : solution)
				solution_prime.push_back(x);

			int v1 = rand() % (problem->n - count_merge);
			int v2 = rand() % (problem->n - count_merge);

			while (v2 == v1)
				v2 = rand() % (problem->n - count_merge);
			// v1 and v2 are different ==> merge partitioning v2 and v1 
			if (v1 > v2)  // make sure that v2 is always bigger than v1 
			{
				int temp = v2;
				v2 = v1;
				v1 = temp;
			}

			for (int i = 0; i < solution_prime.size(); ++i)
			{

				if (solution_prime[i] == v2)
					solution_prime[i] = v1;
				if (solution_prime[i] > v2)
					solution_prime[i]--; // move index of partitioning down 
			}


			double obj_sol_prime = objective(solution_prime, problem->n - count_merge - 1);

			if (obj_sol_prime < best_obj_move)
			{
				best_obj_move = obj_sol_prime;
				move_v1 = v1;
				move_v2 = v2;
			}
			solution_prime.clear();
			move++;
		}

		for (int i = 0; i < solution.size(); ++i)
		{
			if (solution[i] == move_v2)
				solution[i] = move_v1;

			if (solution[i] > move_v2)
				solution[i]--; // move index of partitioning down 
		}

		// next iteration of merging
		count_merge++;

		if (count_merge == problem->n - problem->k)
			break;
	}

	for (auto x : solution)
		cout << x << " ";

	double obj = objective(solution);
	cout << " obj : " << obj << endl;
	// solution rearrangement  neccesary for the CMSA
	rearrange(solution);

	if (problem->alg == KMEANG)
	{
		double end_time = timer.elapsed_time(Timer::VIRTUAL);
		cout << "time: " << (end_time - cur_time) << "\n";
		// store the solution: 
		ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m) +
			"_" + std::to_string(problem->k) + "_" + std::to_string(problem->d) + "_" + std::to_string(problem->alg) + ".out");
		// write  into a file
		//std::cout << std::setprecision(12) << std::fixed;
		myfileOut << "obj: " << obj << endl;
		myfileOut << "value: " << std::fixed << obj << endl;
		myfileOut << "time: " << (end_time - cur_time) << "\n";
		myfileOut.close();
	}
	return solution;
}

// CMSA hybrid algorithm -- helper structures 

struct hash_pair {
	template <class T1, class T2>
	size_t operator()(const pair<T1, T2>& p) const
	{
		auto hash1 = hash<T1>{}(p.first);
		auto hash2 = hash<T2>{}(p.second);

		//if (hash1 != hash2) {
		return hash1 ^ hash2;
		// }

		 // If hash1 == hash2, their XOR is zero.
		 //  return hash1;
	}
};

vector<int> ApplyExactSolver(set<pair<int, int>>& C_prime, double UB = INF)
{
	//cout << "Applying solver..." << endl;
	vector<int> S_sub;
	if (cmsa_milp == 0)
		S_sub = cplex_COAM(C_prime, UB);
	else
		S_sub = cplex_init_faria(C_prime, UB);
	//cout << "Solving finished." << endl;
	return S_sub;
}

void Adapt(unordered_map<pair<int, int>, int, hash_pair>& Age, set<pair<int, int>>& C_prime, vector<int>& S_opt_prime)
{
	for (auto c : C_prime)
	{

		Age[c]++;

	}
	int v_index = 0;
	for (auto j : S_opt_prime)
	{
		pair<int, int> c = make_pair(v_index, j);
		v_index++;
		Age[c] = 0; // promising components  
	}

	for (auto c : C_prime)
	{

		if (Age[c] >= age_max)
		{
			C_prime.erase(c);// cout << "REMOVE COMPONENTS " << c.first << " " << c.second << endl;
			Age.erase(c);
		}

	}
}


// number of appearances of element x in solution
int no_appear(vector<int>& solution, int x)
{

	int counter = 0;
	for (auto el : solution)
		if (x == el)
			counter++;

	return counter;
}


/** LS from GA (Aleksandar's implementation) **/

double move_fit(vector<int>& sol, int i, int p, double** p_sum)
{
	for (int j = 0; j < problem->m; j++) {

		p_sum[sol[i]][j] -= wmat(i, j);
		p_sum[p][j] += wmat(i, j);
	}
	double new_fit = 0;

	for (int z = 1; z < problem->k; z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem->m; q++) {

				double diff = abs(p_sum[z][q] - p_sum[j][q]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;

			}
	}
	for (int j = 0; j < problem->m; j++) {
		p_sum[sol[i]][j] += wmat(i, j);
		p_sum[p][j] -= wmat(i, j);
	}
	return new_fit;
}

double swap_fit(vector<int>& sol, int i, int j, double** p_sum)
{
	for (int s = 0; s < problem->m; s++) {
		p_sum[sol[i]][s] -= wmat(i, s);
		p_sum[sol[j]][s] -= wmat(j, s);
		p_sum[sol[i]][s] += wmat(j, s);
		p_sum[sol[j]][s] += wmat(i, s);
	}
	double new_fit = 0;
	for (int z = 1; z < problem->k; z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem->m; q++) {
				double diff = abs(p_sum[z][q] - p_sum[j][q]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;
			}
	}
	for (int s = 0; s < problem->m; s++) {
		p_sum[sol[i]][s] += wmat(i, s);
		p_sum[sol[j]][s] += wmat(j, s);
		p_sum[sol[i]][s] -= wmat(j, s);
		p_sum[sol[j]][s] -= wmat(i, s);
	}
	return new_fit;
}

double swap_fit3(vector<int>& sol, int i, int j, int q, double** p_sum)
{
	// i gets what j has
	// j gets what q has
	// q gets what i has
	for (int s = 0; s < problem->m; s++) {
		p_sum[sol[i]][s] -= wmat(i, s);
		p_sum[sol[j]][s] -= wmat(j, s);
		p_sum[sol[q]][s] -= wmat(q, s);
		p_sum[sol[j]][s] += wmat(i, s); // i is moved to partition where is currently j
		p_sum[sol[q]][s] += wmat(j, s); // j is moved to partition where is currently q
		p_sum[sol[i]][s] += wmat(q, s); // q is moved to partition where is currently i
	}
	double new_fit = 0;
	for (int z = 1; z < problem->k; z++) {
		for (int j = 0; j < z; j++)
			for (int r = 0; r < problem->m; r++) {
				double diff = abs(p_sum[z][r] - p_sum[j][r]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;
			}
	}
	for (int s = 0; s < problem->m; s++) {
		p_sum[sol[i]][s] += wmat(i, s);
		p_sum[sol[j]][s] += wmat(j, s);
		p_sum[sol[q]][s] += wmat(q, s);
		p_sum[sol[j]][s] -= wmat(i, s);
		p_sum[sol[q]][s] -= wmat(j, s);
		p_sum[sol[i]][s] -= wmat(q, s);
	}
	return new_fit;
}

double move_fit2(vector<int>& sol, int i1, int p1, int i2, int p2, double** p_sum)
{
	if (i1 == i2) {
		cout << "In move_fit2 i and j must be different." << endl;
		exit(1);
	}

	for (int j = 0; j < problem->m; j++) {

		p_sum[sol[i1]][j] -= wmat(i1, j);
		p_sum[p1][j] += wmat(i1, j);
		p_sum[sol[i2]][j] -= wmat(i2, j);
		p_sum[p2][j] += wmat(i2, j);
	}
	double new_fit = 0;

	for (int z = 1; z < problem->k; z++) {
		for (int j = 0; j < z; j++)
			for (int q = 0; q < problem->m; q++) {
				double diff = abs(p_sum[z][q] - p_sum[j][q]); // partial fit. evaluation
				if (diff > new_fit)
					new_fit = diff;

			}
	}
	for (int j = 0; j < problem->m; j++) {
		p_sum[sol[i1]][j] += wmat(i1, j);
		p_sum[p1][j] -= wmat(i1, j);
		p_sum[sol[i2]][j] += wmat(i2, j);
		p_sum[p2][j] -= wmat(i2, j);
	}
	return new_fit;
}

double LSfirst(vector<int>& sol)
{

	int k = problem->k;
	double fit = objective(sol); //problem.fitness(sol, k);
	double** p_sum = new double* [k];

	for (int i = 0; i < k; i++) {
		p_sum[i] = new double[problem->m];
		for (int j = 0; j < problem->m; j++)
			p_sum[i][j] = 0;
	}

	for (int i = 0; i < problem->n; i++)
		for (int j = 0; j < problem->m; j++)
			p_sum[sol[i]][j] += wmat(i, j);

	int impr = 1;
	int nx = problem->n; // n = 1
	while (impr) {
		impr = 0;
		// LS1
		int i_r = rand() % problem->n;
		for (int ix = 0; ix < nx; ix++) {
			int i = (i_r + ix) % problem->n; // to avoid positional bias
			int p = sol[i];
			int oldp = sol[i];
			int p_r = rand() % k;
			for (int px = 0; px < k; px++) {
				//do {
				//	p = rand() % k;
				//} while (p == sol[i]);
				p = (px + p_r) % k;
				if (p == sol[i])
					continue;
				double new_fit = move_fit(sol, i, p, p_sum);

				if (new_fit < fit - 0.1) {
					for (int j = 0; j < problem->m; j++) {

						p_sum[sol[i]][j] -= wmat(i, j);
						p_sum[p][j] += wmat(i, j);
					}

					sol[i] = p;
					fit = new_fit;
					impr = 1;
					//double control_fit = objective(sol, k);
					//if (fabs(control_fit - fit) >1)
					//    cout << "LS1 Incorrect partial fitness function." << to_string(control_fit) << " vs " << to_string(fit) << endl;
					//cout << "LS1 " << to_string(new_fit) << endl;
					break;
				}
			}
			if (impr)
				break;
		}
		if (impr)
			continue;
		// LS2
		int ir = rand() % problem->n;
		for (int ix = 0; ix < nx; ix++) {
			int i = (ir + ix) % problem->n; // to avoid positional bias
			int jr = rand() % problem->n;

			for (int jx = 0; jx < nx; jx++) {
				int j = (jr + jx) % problem->n;
				if (i >= j)
					continue;
				double new_fit = swap_fit(sol, i, j, p_sum);
				if (new_fit < fit - 0.1) {

					for (int s = 0; s < problem->m; s++) {

						p_sum[sol[i]][s] -= wmat(i, s);
						p_sum[sol[j]][s] -= wmat(j, s);
						p_sum[sol[i]][s] += wmat(j, s);
						p_sum[sol[j]][s] += wmat(i, s);
					}
					int pi = sol[i];
					sol[i] = sol[j];
					sol[j] = pi;
					fit = new_fit;
					impr = 1;
					//double control_fit = objective(sol);
					//if (fabs(control_fit - fit) > 1)
					//	cout << "LS2 Incorrect partial fitness function." << to_string(control_fit) << " vs " << to_string(fit) << endl;
					//cout << "LS2 " << to_string(new_fit) << endl;
					break;
				}
			}
			if (impr)
				break;
		}
	}

	for (int i = 0; i < k; i++)
		delete[] p_sum[i];

	return fit;
}

double LSbest(vector<int>& sol)
{
	int n = problem->n;
	int k = problem->k;
	double fit = objective(sol);
	double** p_sum = new double* [k];

	for (int i = 0; i < k; i++) {
		p_sum[i] = new double[problem->m];
		for (int j = 0; j < problem->m; j++)
			p_sum[i][j] = 0;
	}

	for (int i = 0; i < problem->n; i++)
		for (int j = 0; j < problem->m; j++)
			p_sum[sol[i]][j] += wmat(i, j);

	// pre-processing counter of elements appearances 
	unordered_map<int, int> maps_count;
	for (auto& x : sol)
	{
		if (maps_count.find(x) == maps_count.end())
			maps_count.insert({ x, 1 });
		else
			maps_count[x]++;
	}
	// cout << "maps=> 0 " << maps_count[0] << endl;
	int impr = 1;
	while (impr) {
		impr = 0;
		// LS1 best
		int best_i = -1;
		int best_p = -1;
		double best_fit = fit;
		// preprocessing 

		for (int i = 0; i < n; i++) {
			if (maps_count[sol[i]] == 1)
				continue;
			for (int p = 0; p < k; p++) {
				if (p == sol[i])   //count(sol.begin(), sol.end(), sol[i] ) == 1 ) // restrictions to make 1-move
					continue;

				double new_fit = move_fit(sol, i, p, p_sum);

				if (new_fit < best_fit - 0.1) {
					best_i = i;
					best_p = p;
					best_fit = new_fit;
					impr = 1;
				}
			}
		}
		if (impr) {

			// update counter 
			maps_count[sol[best_i]]--;
			maps_count[best_p]++;

			for (int j = 0; j < problem->m; j++) {

				p_sum[sol[best_i]][j] -= wmat(best_i, j);
				p_sum[best_p][j] += wmat(best_i, j);
			}

			sol[best_i] = best_p;
			fit = best_fit;
			impr = 1;
			//double control_fit = objective(sol, k);
			//if (fabs(control_fit - fit) >1)
			//    cout << "LS1 Incorrect partial fitness function." << to_string(control_fit) << " vs " << to_string(fit) << endl;
			//cout << "LS1 " << to_string(fit) << endl;
			continue;
		}
		// LS2 (no need to update mpas_counter, 1-swap) 
		if (ls_number >= 2)
		{
			best_i = -1;
			int best_j = -1;
			best_fit = fit;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < i; j++) {
					if (sol[i] == sol[j])
						continue;
					double new_fit = swap_fit(sol, i, j, p_sum);
					if (new_fit < best_fit - 0.1) {
						best_fit = new_fit;
						best_i = i;
						best_j = j;
						impr = 1;
					}
				}
			}
			if (impr) {
				for (int s = 0; s < problem->m; s++) {

					p_sum[sol[best_i]][s] -= wmat(best_i, s);
					p_sum[sol[best_j]][s] -= wmat(best_j, s);
					p_sum[sol[best_i]][s] += wmat(best_j, s);
					p_sum[sol[best_j]][s] += wmat(best_i, s);
				}
				int pi = sol[best_i];
				sol[best_i] = sol[best_j];
				sol[best_j] = pi;
				fit = best_fit;
				impr = 1;
				//double control_fit = objective(sol);
				//if (fabs(control_fit - fit) > 1)
				//	cout << "LS2a Incorrect partial fitness function." << to_string(control_fit) << " vs " << to_string(fit) << endl;
				//cout << "LS2a " << to_string(fit) << endl;
				continue;
			}
		}
		// LS3 find (i, j, q) and interchange their parts
		if (ls_number >= 3)
		{
			int best_i = -1;
			int best_j = -1;
			int best_q = -1;
			best_fit = fit;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (sol[i] == sol[j])
						continue;
					for (int q = 0; q < n; q++) {
						if (sol[q] == sol[i] || sol[q] == sol[j])
							continue;
						// when k=2 this search does not make sense -- because the next command is never achieved.
						double new_fit = swap_fit3(sol, i, j, q, p_sum);
						if (new_fit < best_fit - 0.1) {
							best_fit = new_fit;
							best_i = i;
							best_j = j;
							best_q = q;
							impr = 1;
						}
					}
				}
			}
			if (impr) {
				// i gets what j has
				// j gets what q has
				// q gets what i has (no changes in maps_count structure)
				for (int s = 0; s < problem->m; s++) {

					p_sum[sol[best_i]][s] -= wmat(best_i, s);
					p_sum[sol[best_j]][s] -= wmat(best_j, s);
					p_sum[sol[best_q]][s] -= wmat(best_q, s);
					p_sum[sol[best_j]][s] += wmat(best_i, s);
					p_sum[sol[best_q]][s] += wmat(best_j, s);
					p_sum[sol[best_i]][s] += wmat(best_q, s);
				}
				int pi = sol[best_i];
				sol[best_i] = sol[best_j];
				sol[best_j] = sol[best_q];
				sol[best_q] = pi;
				fit = best_fit;
				impr = 1;
				double control_fit = objective(sol);
				if (fabs(control_fit - fit) > 1)
					cout << "LS3 Incorrect partial fitness function " << to_string(control_fit) << " vs " << to_string(fit) << endl;
				//cout << "LS3 " << to_string(fit) << endl;
				continue;
			}
		}
	}

	for (int i = 0; i < k; i++)
		delete[] p_sum[i];

	// due to small error accumulation in LS2 we just recalculate objective using standard full objective function 
	fit = objective(sol);

	return fit;
}

/** End of LS from GA **/

// 1-move --- pick an element i and move it to (random) paritioning j 
// restriction: s_j > s_i  ==> and change s_j = s_i
vector<int> LS1(vector<int>& solution)
{

	cout << "Local search... " << endl;
	int iter_no_impr = 0;
	vector<int> best = solution; double best_obj = objective(solution);


	while (true)
	{

		vector<int> s_prime = best;

		int i = rand() % solution.size();
		int j = rand() % solution.size();
		while (s_prime[j] < s_prime[i])
		{
			i = rand() % solution.size();
			j = rand() % solution.size();
		}



		s_prime[j] = s_prime[i];
		double s_prime_obj = objective(s_prime);

		if (s_prime_obj < best_obj) //First improvement strategy
		{
			cout << "improved: " << s_prime_obj << endl;
			best_obj = s_prime_obj;
			best = s_prime;
			iter_no_impr = 0;
		}
		else
			iter_no_impr++;

		if (iter_no_impr > 10000) { // should be parametrized
			cout << "No improvement " << endl;
			break;
		}
	}

	return best;


}

// k--swap based LS 
vector<int> LS2(vector<int>& solution)
{

	cout << "Local search. LS2 .. " << endl;
	int iter_no_impr = 0;
	vector<int> best = solution; double best_obj = objective(solution);


	bool improved = true;

	while (improved)
	{

		improved = false;
		vector<int> s_prime_curr = best;
		for (int i = 0; i < problem->n; ++i)
		{
			for (int j = 0; j < problem->n; ++j)
			{

				vector<int> s_prime = s_prime_curr;
				int k = 0;
				int kx = 1;  // kx-swap -- maximal neighbourhood
				while (k < kx)
				{	// swap s_i and s_j
					if (i != j)
					{
						int temp = s_prime[j];
						s_prime[j] = s_prime[i];
						s_prime[i] = temp;
					}
					++k;

				}

				double s_prime_obj = objective(s_prime);
				cout << "s_prime_obj: " << s_prime_obj << endl;
				if (s_prime_obj < best_obj) //Best improvement strategy
				{
					cout << "solution improved... " << endl;
					best_obj = s_prime_obj;
					best = s_prime;
					iter_no_impr = 0;
					improved = true;
				}
				else
					iter_no_impr++;
			}
		}
	}


	return best;


}

void CMSA()
{

	int cmsa_iter_best = 0; int cmsa_iter = 1;
	unordered_map<pair<int, int>, int, hash_pair> Age;
	set<pair<int, int>> C_prime;

	double end_time = timer.elapsed_time(Timer::VIRTUAL);

	while (end_time - cur_time <= vreme)
	{

		best_rand_sol = INF;

		for (int i = 0; i < n_a; ++i) // generate random solutions:
		{

			vector<int> S;
			if (cmsa_greedy == MDRGHG)
				S = MDRGH(problem->d);
			else
				S = KMeansHeuristic();

			if (objective(S) < best_rand_sol)
				best_rand_sol = objective(S);

			for (int i = 0; i < S.size(); ++i)
			{

				std::pair <int, int> c = make_pair(i, S[i]);

				if (C_prime.find(c) == C_prime.end()) { // not found

					Age[c] = 0;
					C_prime.insert(c);
				}
			}
		}
		vector<int> S_opt_prime;
		S_opt_prime = (ApplyExactSolver(C_prime)); // allow 10 sec for solving a sub-instance ==> should be a parameter
		// call an LS procedure:
		LSbest(S_opt_prime);

		if (!S_opt_prime.empty() and objective(S_opt_prime) < obj_best)
		{
			obj_best = objective(S_opt_prime);
			s_bsf = S_opt_prime;
			cmsa_iter_best = cmsa_iter;
		}

		Adapt(Age, C_prime, S_opt_prime);
		cmsa_iter++;
		end_time = timer.elapsed_time(Timer::VIRTUAL); // update the execution time
		cout << "OBJ BEST: " << obj_best << endl;
	}

	// report stats: 
	cout << "obj : " << obj_best << endl;
	end_time = timer.elapsed_time(Timer::VIRTUAL);
	cout << "time: " << (end_time - cur_time) << endl;
	cout << "iter_best: " << cmsa_iter_best << endl;

	// store the solution: 
	ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m)
		+ "_" + std::to_string(problem->k) + "_" + std::to_string(problem->d)
		+ "_" + std::to_string(problem->alg)
		+ "_" + std::to_string(int(n_a)) + "_" + std::to_string(int(age_max))
		+ "_" + std::to_string(int(cmsa_cplex_time))
		+ "_" + std::to_string(indeks) + ".out");
	// write  into a file
	myfileOut << "value: " << std::fixed << obj_best << endl;
	myfileOut << "time: " << (end_time - cur_time) << endl;
	myfileOut << "iter_best: " << cmsa_iter_best << endl;
	myfileOut.close();


}

/** Adapt-CMSA implementation **/

vector<int> ProbabilisticSolutionConstruction(vector<int>& s_bs, double alpha_bsf)
{

	vector<int> s = s_bs;      //   srand(time(NULL));
	double mutate = 1.0 - alpha_bsf; // probability of mutation
	//cout << "mutate factor " << mutate << endl;
	for (int i = 0; i < s.size(); ++i)
	{

		double r = ((double)rand() / (RAND_MAX)) + 1;
		if (r - 1 < mutate) // do mutation
		{

			/*int oldp = s[i];
			do {
				s[i] = rand() % problem->k;
			} while (s[i] == oldp);*/

			// cout <<to_string(oldp) << " --> "<< to_string(s[i]) << endl;

			int partition = rand() % problem->k; // which partition index to choose

			if (no_appear(s, s[i]) <= 1)
			{
				continue;
			}
			else {
				s[i] = partition;
			}

		}
	}
	// re-arrange parition 
	rearrange(s);
	return s;
}

/** swap-based prob. sol. gen. based on the biased param alpha_bsf  **/
vector<int> ProbabilisticSolutionConstructionSwapBased(vector<int>& s_bs, double alpha_bsf)
{

	vector<int> s = s_bs;      //   srand(time(NULL));
	double mutate = 1.0 - alpha_bsf; // probability of mutation
	//cout << "mutate factor " << mutate << endl;
	for (int i = 0; i < s.size() - 1; ++i)
	{
		for (int j = i + 1; j < s.size(); ++j)
		{

			double r = ((double)rand() / (RAND_MAX)) + 1;
			if (r - 1 < mutate) // do mutation
			{
				int part1 = s[i];
				int part2 = s[j];
				s[i] = part2;
				s[j] = part1;
			}
		}
	}
	// re-arrange parition (if nencessary)
	rearrange(s);
	return s;
}


/** solve subinstance partially ==> inclusivelly  **/
// input: set of  components: C_prime 

struct ComparePairs {

	ComparePairs(unordered_map<pair<int, int>, int, hash_pair>& Appearances) { this->Appearances = Appearances; }
	bool operator () (pair<int, int>& i, pair<int, int>& j)
	{

		if (Appearances[i] > Appearances[j])
			return 1;
		if (Appearances[i] <= Appearances[j])
			return 0;

	}

	unordered_map<pair<int, int>, int, hash_pair>  Appearances;

};

vector<int> formSubinstanceDeep(std::vector<vector<pair<int, int>>>& C_prime)
{

	vector<set<int>> Account;
	// init
	for (int i = 0; i < problem->n; ++i)
	{
		set<int> s_i;
		Account.push_back(s_i);
	}
	//fill the vector
	for (int i = 0; i < C_prime.size(); ++i)
		for (auto& x : C_prime[i])
			Account[x.first].insert(x.second);

	/*for(int i = 0; i < problem->n; ++i)
		cout << Account[i].size() << endl;*/

	vector<int> sort_array;
	for (int i = 0; i < problem->n; ++i)
		sort_array.push_back(i);

	//find the sorted order acc. to Account 
	for (int i = 0; i < problem->n - 1; ++i)
		for (int j = i + 1; j < problem->n; ++j)
		{
			if (Account[sort_array[i]].size() < Account[sort_array[j]].size())
			{
				int temp = sort_array[i];
				sort_array[i] = sort_array[j];
				sort_array[j] = temp;
			}
		}

	// generate a series of subinstances and solve them 
	vector<int> best_sol_iter;  double sol_in_iter = INF;  set<pair<int, int>> C_to_j;  // C_j: subinstance per iteration (increassing)

	int start = problem->k / 2; //should be parametrized

	//conversion sets to vectors
	vector<vector<int>> Account_vector;
	for (int i = 0; i < problem->n; ++i)
	{
		std::vector<int> v(Account[i].begin(), Account[i].end());
		Account_vector.push_back(v);
	}

	for (int i = 0; i < problem->n; ++i) // first start-elements in each 
	{
		for (int k = 0; k < start; ++k)
		{

			if (Account_vector[sort_array[i]].size() > k)
			{
				pair<int, int> px = std::make_pair(sort_array[i], Account_vector[sort_array[i]][k]);
				C_to_j.insert(px);
			}
		}
	}
	// solve the initial subinstance
	vector<int> S_bsf_iter;
	S_bsf_iter = ApplyExactSolver(C_to_j, sol_in_iter);

	if (!S_bsf_iter.empty())
	{

		sol_in_iter = objective(S_bsf_iter); // update sol_in_iter ==> always the best along the iterations
		best_sol_iter = S_bsf_iter;

	}

	for (int l = start; l < problem->k; ++l) // k ==> k + 1
	{
		cout << "l-th step =====================================> " << l << endl;
		for (int v = 0; v < problem->n; ++v) // first start-elements in each 
		{
			//cout << "sort_array[ v ]" << sort_array[ v ] << endl;
			if (Account_vector[sort_array[v]].size() > l)
			{
				pair<int, int> px = std::make_pair(sort_array[v], Account_vector[sort_array[v]][l]);
				C_to_j.insert(px);
			}
		}
		// solve adapted subinstance 
		vector<int> S_opt_prime;
		S_opt_prime = ApplyExactSolver(C_to_j, sol_in_iter);

		// update sols
		if (!S_opt_prime.empty()) // keep best solution if a local solution is found
		{
			// if within 0.9 best solution ==> improve by the LS
			double  s_opt_prime_obj = objective(S_opt_prime);
			/*if(s_opt_prime_obj >= 0.6 * objective(s_bsf))
			{
			   cout << "Apply LS procedure ======================================================>  " << endl;
			   s_opt_prime_obj = LS(S_opt_prime); // possibly improve
			}*/

			if (s_opt_prime_obj < sol_in_iter)
			{
				best_sol_iter.clear();

				for (auto x : S_opt_prime)
					best_sol_iter.push_back(x);

				sol_in_iter = objective(S_opt_prime); // update sol_in_iter ==> always the best along the iterations
			}
		}

		double end_time = timer.elapsed_time(Timer::VIRTUAL);

		if (end_time - cur_time >= vreme) // time break after time limit has exceeded;
			break;
	}
	return best_sol_iter;

}


vector<int> formSubinstance(std::vector<vector<pair<int, int>>>& C_prime)
{
	//cout << "Form subinstance based on a series of subinstances and its UBs... " << C_prime.size() << endl;
	unordered_map<pair<int, int>, int, hash_pair> Appearances;


	for (int i = 0; i < C_prime.size(); ++i)
	{
		for (std::vector<pair<int, int>>::iterator it = C_prime[i].begin(); it != C_prime[i].end(); ++it)
		{

			if (Appearances.find(*it) == Appearances.end()) // still not exists in Appeances
				Appearances[*it] = 1;
			else
				Appearances[*it]++;
		}
	}
	// sort pairs (components) 

	for (int i = 0; i < C_prime.size(); ++i)

		sort(C_prime[i].begin(), C_prime[i].end(), ComparePairs(Appearances));



	vector<int> best_sol_iter;  double sol_in_iter = INF;  set<pair<int, int>> C_to_j;

	int addFirst = problem->n / 2; // this has to be a parameter

	for (int j = 0; j < addFirst; ++j)
	{
		for (int i = 0; i < C_prime.size(); ++i)
			C_to_j.insert(C_prime[i][j]);
	}
	vector<int> S_bsf_iter;
	S_bsf_iter = ApplyExactSolver(C_to_j, sol_in_iter); // allow 10 sec for solving a sub-instance ==> should be a parameter

	if (!S_bsf_iter.empty())
	{

		sol_in_iter = objective(S_bsf_iter); // update sol_in_iter ==> always the best along the iterations
		best_sol_iter = S_bsf_iter;

	}
	// add one-by-one
	for (int j = addFirst; j < C_prime[0].size(); ++j) // i = addFrist, ..., n (partially)
	{
		//cout << "added most occured j-th column " << j << endl; 
		vector<int> S_opt_prime;

		for (int i = 0; i < C_prime.size(); ++i)
			C_to_j.insert(C_prime[i][j]);

		S_opt_prime = ApplyExactSolver(C_to_j, sol_in_iter); // allow 10 sec for solving a sub-instance ==> should be a parameter

		if (!S_opt_prime.empty() and objective(S_opt_prime) < sol_in_iter) // keep best solution
		{
			best_sol_iter.clear();

			for (auto x : S_opt_prime)
				best_sol_iter.push_back(x);

			sol_in_iter = objective(S_opt_prime); // update sol_in_iter ==> always the best along the iterations
		}

		double end_time = timer.elapsed_time(Timer::VIRTUAL);

		if (end_time - cur_time >= vreme) // time break after time limit has exceeded;
			break;
	}

	return best_sol_iter;
}


void Adapted_CMSA()
{
	double eps = 1;
	int cmsa_iter_best = 0; int cmsa_iter = 0; int ls_iter_impr = 0;
	unordered_map<pair<int, int>, int, hash_pair> Age;
	set<pair<int, int>> C_prime;

	int n_a_init = n_a; //n_a initialized from a command line  (n_a = 1, in the case of the original version of AdaptCMSA)
	double alpha_bsf = alpha_UB;

	/** start the procedure **/
	if (cmsa_greedy == MDRGHG)
		s_bsf = MDRGH(problem->d);
	else
	if (cmsa_greedy == KMEANG)
	    s_bsf = KMeansHeuristic();
	else
	    s_bsf = RandomizedGenerator();


	obj_best = objective(s_bsf);

	int  ix = 0;
	vector<vector<pair<int, int>>> C_prime_vector;   // for Deep-CMSA
	vector<pair<int, int>> C_0;
	for (auto x : s_bsf)
	{
		std::pair <int, int> c = make_pair(ix, x);
		C_prime.insert(c);

		if (problem->alg == DEEP_CMSAH) //Deep-CMSA
			C_0.push_back(c);

		ix++;
	}
	if (problem->alg == DEEP_CMSAH) //Deep-CMSA
		C_prime_vector.push_back(C_0);// first vector to index 0

	double end_time = timer.elapsed_time(Timer::VIRTUAL);
	double S_opt_prime_obj = obj_best;
	// start with the main iterations 
	while (end_time - cur_time <= vreme)
	{
		/*cout << to_string(cmsa_iter) << ".\tn_a_init=" << to_string(n_a_init) << "\talpha_bsf=" << to_string(alpha_bsf)
			<< "\tbest=" << to_string(obj_best) << "\tnew=" << to_string(S_opt_prime_obj) << endl;*/

		for (int i = 0; i < n_a_init; ++i) // generate random solutions:
		{

			vector<pair<int, int>> C_i;

			vector<int> S;
			//float r = ((double) rand() / (RAND_MAX)) + 1;  

			S = ProbabilisticSolutionConstruction(s_bsf, alpha_bsf);

			double S_obj = LSbest(S);
			rearrange(S);
			if (S_obj < obj_best) // update best sol
			{
				s_bsf.clear();
				for (auto x : S) // update s_bsf
					s_bsf.push_back(x);

				obj_best = S_obj;
			}

			// Merge -- a standard way
			for (int i = 0; i < S.size(); ++i)
			{
				std::pair <int, int> c = make_pair(i, S[i]);
				if (C_prime.find(c) == C_prime.end())// component c not found
				{
					C_prime.insert(c); //  merge in the standard CMSA  
				}
				if (problem->alg == DEEP_CMSAH)
					C_i.push_back(c);
			}
			if (problem->alg == DEEP_CMSAH) //Deep-CMSA
				C_prime_vector.push_back(C_i);

			S.clear(); C_i.clear();
		}
		// solve subinstance:
		double start_time_subinstance = timer.elapsed_time(Timer::VIRTUAL);

		vector<int> S_opt_prime;
		if (problem->alg == DEEP_CMSAH) //Deep-CMSA
			S_opt_prime = formSubinstanceDeep(C_prime_vector);
		//S_opt_prime = formSubinstance( C_prime_vector ); 
		else
			S_opt_prime = ApplyExactSolver(C_prime);

		double end_time_subinstance = timer.elapsed_time(Timer::VIRTUAL);

		// Start of LS application: 
		S_opt_prime_obj = INF;
		if (!S_opt_prime.empty()) // if not the empty solution
		{
			S_opt_prime_obj = objective(S_opt_prime);
			//if (true  || S_opt_prime_obj >= 0.8 * obj_best)
			//{
			ls_iter_impr++;
			S_opt_prime_obj = LSbest(S_opt_prime); // possibly improve
			rearrange(S_opt_prime);
			//double stest = LS(S_opt_prime);
			//if (fabs(stest - S_opt_prime_obj)>0.01)
			//    cout << "GRESKA" << endl;
		//}
		}

		// End of LS application
		double t_solve = end_time_subinstance - start_time_subinstance;

		if (t_solve < t_prop * cmsa_cplex_time and alpha_bsf > alpha_LB)
			alpha_bsf -= alpha_red;

		if (!S_opt_prime.empty() and S_opt_prime_obj < obj_best)
		{
			s_bsf.clear();
			for (auto x : S_opt_prime) // update s_bsf
				s_bsf.push_back(x);

			obj_best = S_opt_prime_obj;
			n_a_init = n_a;
			alpha_bsf = alpha_UB; // aca dodao
			cmsa_iter_best = cmsa_iter;
		}
		else {
			/*
			if (t_solve < t_prop * cmsa_cplex_time) {
				// there was enough time for cplex to deal with model
				// so we adjust gradually alpha_bs to increase model search space
				if(alpha_bsf - alpha_red >= alpha_LB)
					alpha_bsf -= alpha_red;
				else {
					// or if alpha_bs is already difficult, we jump to the higher neighborhood with the easiest alpha (this is gradual change of pair (n_a_init, alpha_bsf)
					n_a_init++;
					alpha_bsf = alpha_UB;
				}
			}
			else {
				// opposite to previous, gradual decrease of model search space through dynamics of (n_a_init, alpha_bsf)
				if (alpha_bsf < alpha_UB)
					alpha_bsf = min(alpha_bsf + alpha_red / 10, alpha_UB);
				else if (n_a_init > n_a){
					n_a_init--;
					alpha_bsf = alpha_LB;
				}
			}*/


			if (!S_opt_prime.empty() and S_opt_prime_obj > obj_best + eps) // if empty ==> the same solution 
			{

				if (n_a_init == n_a)
					alpha_bsf = min(alpha_bsf + alpha_red / 10, alpha_UB);
				else {
					n_a_init = n_a;
					alpha_bsf = alpha_UB; // aca dodao
				}
			}
			else   // no improvement => local optima ==> diverse the search for more 
				n_a_init++;

		}
		C_prime.clear();

		if (problem->alg == DEEP_CMSAH)
		{
			for (auto x : C_prime_vector)
				x.clear();
		}

		ix = 0;

		for (auto x : s_bsf)
		{
			std::pair <int, int> c = make_pair(ix, x);
			if (problem->alg == DEEP_CMSAH)
				C_0.push_back(c);
			else
				C_prime.insert(c);

			ix++;
		}
		if (problem->alg == DEEP_CMSAH)
			C_prime_vector.push_back(C_0);

		cmsa_iter++;
		end_time = timer.elapsed_time(Timer::VIRTUAL);

	}

	// report stats:        
	bool valid = validityCheck(s_bsf);
	cout << "value: " << std::fixed << obj_best << endl;
	cout << "time: " << (end_time - cur_time) << endl;
	cout << "CMSA iter best: " << cmsa_iter_best << endl;
	cout << "CMSA iter overall: " << cmsa_iter << endl;
	cout << "LS iter improved: " << ls_iter_impr << endl;
	cout << "validity: " << valid << endl;
	cout << "solution: ";
	for (auto x : s_bsf)
		cout << x << " ";

	// store the solution: 
	ofstream myfileOut(outPath + std::to_string(problem->n) + "_" + std::to_string(problem->m)
		+ "_" + std::to_string(problem->k) + "_" + std::to_string(problem->d)
		+ "_" + std::to_string(problem->alg)
		+ "_" + std::to_string(int(n_a)) + "_" + std::to_string(int(age_max))
		+ "_" + std::to_string(int(cmsa_cplex_time))
		+ "_" + std::to_string(indeks) + ".out");
	// write  into a file
	myfileOut << "value: " << std::fixed << obj_best << endl;
	myfileOut << "time: " << (end_time - cur_time) << endl;
	myfileOut << "CMSA iter best: " << cmsa_iter_best << endl;
	myfileOut << "CMSA iter overall: " << cmsa_iter << endl;
	myfileOut << "LS iter improved: " << ls_iter_impr << endl;
	myfileOut << "validity: " << valid << endl;
	myfileOut << "solution: ";
	for (auto x : s_bsf)
		myfileOut << x << " ";

	myfileOut.close();
}

pair<double, vector<int>> random_solution() {
	// generating random solution
	pair<double, vector<int>> sol;
	sol.second = vector<int>(problem->n);
	for (int j = 0; j < problem->n; j++)
		sol.second[j] = rand() % problem->k;
	rearrange(sol.second);
	sol.first = objective(sol.second);
	return sol;
}

pair<double, vector<int>> select_tournament(vector<pair<double, vector<int>>> pop, int t_size) {

#ifdef _WIN32	
	double best = DBL_MAX;
#else
	double best = std::numeric_limits<double>::max();

#endif
	double best_i = -1;
	for (int i = 0; i < t_size; i++) {
		int ri = rand() % pop.size();
		if (pop[ri].first < best) {
			best = pop[ri].first;
			best_i = ri;
		}
	}
	return pop[best_i];
}
 
 
 
 
/** end of Adapt-CMSA imeplementation **/

/** terminal params **/
void read_parameters(int argc, char** argv) {

	if (argc < 5)
	{
		printf("Not enough arguments \n cplex_mdmwnpp <input_file> <n> <m> <k> <d> <alg> <move> <time> \n \n If <time> == -1, it is unbounded! \n \n ");
		exit(0);
	}
	// problem allocation 
	problem = (problem_struct*)malloc(sizeof(problem_struct));

	int iarg = 1;
	//cout << "Reading params..." << endl;
	while (iarg < argc) {

		if (strcmp(argv[iarg], "-f") == 0)          sprintf(imeul, "%s", argv[++iarg]);
		else if (strcmp(argv[iarg], "-n") == 0)      problem->n = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-m") == 0)      problem->m = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-k") == 0)      problem->k = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-d") == 0)      problem->d = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-alg") == 0)    problem->alg = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-iter_move") == 0)   iter_move = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-t") == 0)      vreme = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-n_a") == 0)  n_a = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-age_max") == 0)  age_max = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-cmsa_cplex_time") == 0)  cmsa_cplex_time = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-cmsa_milp") == 0)  cmsa_milp = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-cmsa_greedy") == 0)  cmsa_greedy = atoi(argv[++iarg]) + 2;
		else if (strcmp(argv[iarg], "-alphaLB") == 0)  alpha_LB = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-alphaUB") == 0)  alpha_UB = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-alpha_red") == 0)  alpha_red = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-t_prop") == 0)   t_prop = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-seed") == 0)   seed = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-index") == 0)  indeks = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-out") == 0)   outPath = argv[++iarg];//sprintf(outPath,"%s",argv[++iarg]);
		else if (strcmp(argv[iarg], "-move_prob") == 0) move_prob = atof(argv[++iarg]);
		else if (strcmp(argv[iarg], "-min_p") == 0) min_p = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-max_p") == 0) max_p = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg], "-ls_number") == 0) ls_number = atoi(argv[++iarg]);

		else ++iarg;
	}
	//cout << "n=" << to_string(problem->n) << " m=" << to_string(problem->m) << " k=" << to_string(problem->k) << endl;

	if (problem->m < problem->d)
		problem->d = problem->m;

	if (problem->d == 0 or problem->d == NULL) // not set 
		problem->d = problem->m;

}

int main(int argc, char** argv)
{
	srand(seed); //time(NULL));
	std::cout << std::setprecision(4) << std::fixed;
	cur_time = timer.elapsed_time(Timer::VIRTUAL);
	read_parameters(argc, argv);
	//cout << "Read from file ..." << endl;
	ulazpod();

	switch (problem->alg)
	{
		/** two relevant models for MDMWNPP from literature **/
	case 0: { cout << "COAM model Nikolic et. al. " << endl;  set<pair<int, int>> C_prime; cplex_COAM(C_prime);                     break;                 }
	case 1: { cout << "Faria et al. model (2020)  " << endl;  set<pair<int, int>> C_prime; cplex_init_faria(C_prime);               break;        	    }
	case 2: { cout << "Greedy: MDRGH " << endl;  MDRGH(problem->d);              break;        	    }
	case 3: { cout << "Greedy: KMeans-based Heuristic " << endl;  KMeansHeuristic(iter_move);     break;        	    }
	case 4: { cout << "CMSA... " << endl;  CMSA();                          break;                 }
	case 5: { /*cout << "Adapted-CMSA " << endl;*/  Adapted_CMSA();                   break;                 }
	case 6: { cout << "Deep-CMSA " << endl;  Adapted_CMSA();                   break;                 }
 
	default: { cout << "Kojic MIP model (2010) " << endl;  cplex_init_kojic();               break;                 }

	}
	//Example of a call: ./mdmwnpp -f mdtwnpp_500_20a.txt -n 50 -m 4 -k 3 -alg 5 -cmsa_cplex_time 3  -cmsa_greedy 2 -cmsa_milp 0 -n_a 1 -alphaLB 0.5 -alphaUB 0.9 -alpha_red 0.05 -t_prop 0.4

	return 0;
} /* main */
