#include <iostream>
#include "Problem.h"
#include "VNS.h"
#include <random>
#include "GA.h"
#include <cstring>

string instance="-1", out="-1", alg="";
int n=-1, m=-1, k=-1,max_seconds=-1, r_max=-1, pop_size=-1, seed=-1, ls3_active=-1;
double crossover_rate=-1, mutation_rate=-1;
Problem problem;

void error_ga(string msg) {
    cout << msg << endl;
    cout << "Expected arguments: ga -f [instance file] -n [n?] -m [m?] -k [k?] -t [max time (s)] -pop_size [population size] -out [output directory] -cross_rate [crossover_rate] -mut_rate [mutation_rate] -ls3_active [ls3 active] -seed [random seed?]" << endl;
    exit(1);
}

void error_vns(string msg) {
    cout << msg << endl;
    cout << "Expected arguments: vns -f [instance file] -n [n?] -m [m?] -k [k?] -t [max time (s)] -r_max [max neighborhood (1-3)] -out [output directory] -seed [random seed?]" << endl;
    exit(1);
}

void error(string msg) {
    cout << msg << endl;
}

/** terminal params **/
void read_parameters(int argc, char** argv) {
    //string alg = argv[1];

    int iarg = 1;
    //cout << "Reading common params..." << endl;
    while (iarg < argc) {
        if (strcmp(argv[iarg], "--alg") == 0)        alg = argv[++iarg];
        else if (strcmp(argv[iarg], "--f") == 0)      instance = argv[++iarg];
        else if (strcmp(argv[iarg], "--n") == 0)      n = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg], "--m") == 0)      m = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg], "--k") == 0)      k = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg], "--t") == 0)      max_seconds = atof(argv[++iarg]);
        else if (strcmp(argv[iarg], "--seed") == 0)   seed = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg], "--out") == 0)    out = argv[++iarg];
        else ++iarg;
    }
    // excplicitly set seed
    if (seed != -1) {
        //cout << "Setting explicit random seed " << to_string(seed) << endl;
        srand(seed);
    }
    if (instance=="-1")
        error("Parameter -f is required.");
    /*if (out=="-1")
        error("Parameter -out is required.");*/

    if (n != -1 && m != -1 && k!=-1)
        problem.load_from_file(instance, n, m, k);
    else
        problem.load_from_file(instance);

    if (alg == "ga") {
        iarg = 1;
        //cout << "Reading ga params..." << endl;
        while (iarg < argc) {

            if (strcmp(argv[iarg], "--pop_size") == 0)    pop_size = (int)(atof(argv[++iarg]) * problem.get_n());
            else if (strcmp(argv[iarg], "--cross_rate") == 0)      crossover_rate = atof(argv[++iarg]);
            else if (strcmp(argv[iarg], "--mut_rate") == 0)      mutation_rate = atof(argv[++iarg]);
            else if (strcmp(argv[iarg], "--ls3_active") == 0)      ls3_active = atoi(argv[++iarg]);
            else ++iarg;
        }
        if(pop_size==-1)
            error_ga("Parameter -pop_size is required.");
        if (crossover_rate == -1)
            error_ga("Parameter -cross_rate is required.");
        if (mutation_rate == -1)
            error_ga("Parameter -mut_rate is required.");
        if (ls3_active == -1)
            error_ga("Parameter -ls3_active is required.");
    }
    else if (alg == "vns") {
        iarg = 1;
        //cout << "Reading ga params..." << endl;
        while (iarg < argc) {

            if (strcmp(argv[iarg], "--r_max") == 0)    r_max = atoi(argv[++iarg]);
            else ++iarg;
        }
        if (r_max == -1)
            error_ga("Parameter -r_max is required.");
    }
    else {
        cout << "Algorithm " << alg << " is not supported." << endl;
        exit(1);
    }
}

int main(int argc, char* argv[])
{
    read_parameters(argc, argv);
    //string alg = argv[1];

    if (alg == "ga") {
        string out_file = out + to_string(problem.get_n()) + "_" + to_string(problem.get_m()) + "_" + to_string(problem.get_k()) + "_" + to_string(seed) + "_" + to_string(crossover_rate)
            + "_" + to_string(mutation_rate) + "_" + to_string(ls3_active) + "_" + alg + ".out";
        GA ga(problem, max_seconds, pop_size, crossover_rate, mutation_rate, ls3_active);
        ga.run();
        ga.output(out_file);
    }
    else if (alg == "vns") {
        string out_file = out + to_string(problem.get_n()) + "_" + to_string(problem.get_m()) + "_" + to_string(problem.get_k()) + "_" + to_string(seed) + "_" + alg + ".out";
        VNS vns(problem, max_seconds, r_max);
        vns.run();
        vns.output(out_file);
    }
    else
        cout << "Not allowed algorithm." << endl;
}
