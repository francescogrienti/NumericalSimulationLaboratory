/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <tuple>
#include "random.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]) {

    //Variables declaration
    Random rnd;
    int M = 100000; //Number of total throws
    int N = 100; //Number of blocks
    double mu = 0.5; //Expected mean
    int expecValue = M / N;
    tuple<vector<double>, vector<double>> averages;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives;
    vector<double> chi2(N, 0.);
    vector<int> count(N, 0);
    int seed[4];
    int p1, p2;

    //Functions
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    //Evaluation of the mean
    averages = mean(M, N, rnd);
    cumulatives = cumulativeAverage(get<0>(averages), get<1>(averages));
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_1.dat");

    //Evaluation of the standard deviation
    averages = mean(M, N, rnd, mu);
    cumulatives = cumulativeAverage(get<0>(averages), get<1>(averages));
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_2.dat");

    //Chi-2
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < M; k++) {
            double r = rnd.Rannyu();
            count[(int) (r * N)] += 1;
        }
        for (int j = 0; j < N; j++) {
            chi2[j] += double(pow(count[j] - expecValue, 2)) / (double) (expecValue);
            count[j] = 0;
        }
    }

    writeOnFile(chi2, "results_3.dat");

    rnd.SaveSeed();
    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
