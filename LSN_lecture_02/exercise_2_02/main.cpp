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
    vector<double> average(N, 0.);
    vector<vector<double>> sum_d(N, vector<double>(L, 0.));
    vector<vector<double>> sum_c(N, vector<double>(L, 0.));
    vector<vector<double>> ave(N, vector<double>(L, 0.));
    vector<vector<double>> ave2(N, vector<double>(L, 0.));
    vector<double> error(N, 0.);
    vector<int> count_d(n, 0);
    vector<int> count_c(n, 0);
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
                sum_d[k][i] += sqrt((pow(count_d[0], 2) + pow(count_d[1], 2) + pow(count_d[2], 2)));
            }
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < L; i++) {
            ave[k][i] = (sum_d[k][i] / L);
            ave2[k][i] = pow(ave[k][i], 2);
        }
    }
    for (int i = 0; i < L; i++) {
        double sum1 = 0.;
        double sum2 = 0.;
        for (int k = 0; k < N; k++) {
            sum1 += ave2[k][i];
            sum2 += ave[k][i];
        }
        average[i] = sum2 / N;
        error[i] = sqrt(((sum1 / N) - pow(sum2 / N, 2)) / N);
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << average[i] << " " << error[i] << "\t" << endl;
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
                sum_c[k][i] += sqrt((pow(count_c[0], 2) + pow(count_c[1], 2) + pow(count_c[2], 2)));
            }
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < L; i++) {
            ave[k][i] = (sum_c[k][i] / L);
            ave2[k][i] = pow(ave[k][i], 2);
        }
    }
    for (int i = 0; i < L; i++) {
        double sum1 = 0.;
        double sum2 = 0.;
        for (int k = 0; k < N; k++) {
            sum1 += ave2[k][i];
            sum2 += ave[k][i];
        }
        average[i] = sum2 / N;
        error[i] = sqrt(((sum1 / N) - pow(sum2 / N, 2)) / N);
    }

    ofstream WriteResults2;
    WriteResults2.open("results_2.dat");
    if (WriteResults2.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << average[i] << " " << error[i] << "\t" << endl;
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
