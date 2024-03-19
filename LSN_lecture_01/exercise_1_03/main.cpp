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
#include "random.h"
#include <vector>
#include <tuple>

double error(std::vector<double> av, std::vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    const int M = 100000; //Throws
    const int N = 100; //Blocks
    const int L = M / N; //Throws per block
    const double l = 1.; //Length of the needle
    const double d = 2.; //Spacing among lines
    int N_hit = 0;
    vector<double> ave(N, 0.);
    vector<double> ave2(N, 0.);
    vector<double> sum_prog(N, 0.);
    vector<double> sum2_prog(N, 0.);
    vector<double> err_prog(N, 0.);
    //tuple<vector<double>, vector<double >> averages;
    //tuple<vector<double>, vector<double>, vector<double >> cumulatives;
    int seed[4];
    int p1, p2;
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

    //Block of code for the evaluation of the mean
    for (int i = 0; i < N; i++) {
        N_hit = 0;
        for (int j = 0; j < L; j++) {
            double x1 = rnd.Rannyu(0, d);
            double x2 = 0.;
            double x = 0;
            double y = 0.;
            double angle = 0.;
            do {
                y = rnd.Rannyu(-1,1);
                x = rnd.Rannyu(-1,1);
            } while (pow(x, 2) + pow(y, 2) > 1.);
            angle = atan(y/x);
            x2 = x1 + l * cos(2 * angle);
            if (x2 < 0 || x2 > d) {
                N_hit += 1;
            }
            ave[i] = (2.0 * l * L) / (N_hit * d); //Store average values for each block
            ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
        }
    }

    for (int k = 0; k < N; k++) {
        sum_prog[k] = 0.;
        sum2_prog[k] = 0.;
        for (int l = 0; l < k + 1; l++) {
            sum_prog[k] += ave[l];
            sum2_prog[k] += ave2[l];
        }
        sum_prog[k] /= (k + 1); //Cumulative average
        sum2_prog[k] /= (k + 1); //Cumulative square average
        err_prog[k] = error(sum_prog, sum2_prog, k); //Statistical uncertainty
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << sum_prog[i] << " " << err_prog[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults1.close();
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
