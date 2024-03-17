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

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    int M = 10000; //Number of throws
    int n = 3; //Possible directions (both positive and negative)
    int N = 100; //Number of blocks
    int L = M / N; //Number of throws per block;
    double a = 1.; //Lattice step
    double *sum_d = new double[N]{0.};
    double *sum_c = new double[N]{0.};
    double *sum_err_d = new double[N]{0.};
    double *sum_err_c = new double[N]{0.};
    vector<vector<double>> err_d(M, vector<double>(N, 0.));
    vector<vector<double>> err_c(M, vector<double>(N, 0.));
    double *error_d = new double[N]{0.};
    double *error_c = new double[N]{0.};
    int *count_d = new int[n]{0};
    double *count_c = new double[n]{0};
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

    //Discrete case
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            count_d[0] = 0;
            count_d[1] = 0;
            count_d[2] = 0;
            for (int i = 0; i < L; i++) {
                double r = rnd.Rannyu(-1, 1);
                if (r < 0) {
                    count_d[(int) ((-1) * r * n)] -= a;
                } else if (r > 0) {
                    count_d[(int) (r * n)] += a;
                }
                sum_d[i] += (pow(count_d[0], 2) + pow(count_d[1], 2) + pow(count_d[2], 2));
                err_d[j + k * N][i] = (pow(count_d[0], 2) + pow(count_d[1], 2) + pow(count_d[2], 2));
            }
        }
    }
    for (int k = 0; k < N; k++) {
        sum_d[k] = sqrt(sum_d[k] / M);
    }

    for (int i = 0; i < N; i++) {
        double somma = 0.;
        for (int k = 0; k < M; k++) {
            somma += pow((sqrt(err_d[k][i]) - sum_d[i]), 2);
        }
        error_d[i] = sqrt(somma / M);
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << sum_d[i] << " " << error_d[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults1.close();

    //Continuum case
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            count_c[0] = 0.;
            count_c[1] = 0.;
            count_c[2] = 0.;
            for (int i = 0; i < L; i++) {
                double theta = rnd.Rannyu(0, M_PI);
                double phi = rnd.Rannyu(0, 2 * M_PI);
                count_c[0] += sin(theta) * sin(phi);
                count_c[1] += sin(theta) * cos(phi);
                count_c[2] += cos(theta);
                sum_c[i] += (pow(count_c[0], 2) + pow(count_c[1], 2) + pow(count_c[2], 2));
                err_c[j + k * N][i] = pow(count_c[0], 2) + pow(count_c[1], 2) + pow(count_c[2], 2);
            }
        }
    }

    for (int k = 0; k < N; k++) {
        sum_c[k] = sqrt(sum_c[k] / M);
    }

    for (int i = 0; i < N; i++) {
        double somma = 0.;
        for (int k = 0; k < M; k++) {
            somma += pow((err_c[k][i] - sum_c[i]), 2);
        }
        error_c[i] = sqrt(somma / M);
    }

    ofstream WriteResults2;
    WriteResults2.open("results_2.dat");
    if (WriteResults2.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults2 << sum_c[i] << " " << error_c[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults2.close();

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
