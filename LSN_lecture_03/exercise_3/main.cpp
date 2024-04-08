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
#include "random.h"
#include "functions.h"
#include "cmath"
#include "algorithm"
#include <vector>

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
    int M = 100000;
    int N = 100;
    int L = M / N;
    int T = 1;
    int strike_price = 100;
    double r = 0.1;
    double sigma = 0.25;
    int S_0 = 100;
    double t_i = (double) T / N;
    vector<double> S(N + 1, 0.);
    vector<double> ave_C_i(N, 0.);
    vector<double> ave2_C_i(N, 0.);
    vector<double> ave_P_i(N, 0.);
    vector<double> ave2_P_i(N, 0.);
    vector<double> sumC_prog(N, 0.);
    vector<double> sum2C_prog(N, 0.);
    vector<double> errC_prog(N, 0.);
    vector<double> sumP_prog(N, 0.);
    vector<double> sum2P_prog(N, 0.);
    vector<double> errP_prog(N, 0.);
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //Direct
    // aggiungere il codice già fatto


    //Discretized
    for (int i = 0; i < N; i++) {
        double sum_C = 0.;
        double sum_P = 0.;
        for (int k = 0; k < L; k++) {
            double z_i = rnd.Gauss(0., 1.);
            S[0] = S_0;
            S[1] = S_0 * exp(((r - pow(sigma, 2) / 2)) * t_i + sigma * z_i * sqrt(t_i));
            for (int j = 1; j < N; j++) {
                z_i = rnd.Gauss(0., 1.);
                S[j + 1] = S[j] * exp(((r - pow(sigma, 2) / 2)) * t_i + sigma * z_i * sqrt(t_i));
            }
            sum_C += exp((-1) * r * T) * max(0., S[N] - strike_price);
            sum_P += exp((-1) * r * T) * max(0., strike_price - S[N]);
        }
        ave_C_i[i] = sum_C / L;
        ave_P_i[i] = sum_P / L;
        ave2_C_i[i] = pow(ave_C_i[i], 2);
        ave2_P_i[i] = pow(ave_P_i[i], 2);
    }

    for (int k = 0; k < N; k++) {
        sumC_prog[k] = 0.;
        sum2C_prog[k] = 0.;
        sumP_prog[k] = 0.;
        sum2P_prog[k] = 0.;
        for (int l = 0; l < k + 1; l++) {
            sumC_prog[k] += ave_C_i[l];
            sum2C_prog[k] += ave2_C_i[l];
            sumP_prog[k] += ave_P_i[l];
            sum2P_prog[k] += ave2_P_i[l];
        }
        sumC_prog[k] /= (k + 1); //Cumulative average
        sum2C_prog[k] /= (k + 1); //Cumulative square average
        errC_prog[k] = error(sumC_prog, sum2C_prog, k); //Statistical uncertainty
        sumP_prog[k] /= (k + 1); //Cumulative average
        sum2P_prog[k] /= (k + 1); //Cumulative square average
        errP_prog[k] = error(sumP_prog, sum2P_prog, k); //Statistical uncertainty
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << sumC_prog[i] << " " << errC_prog[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;

    WriteResults1.close();

    ofstream WriteResults2;
    WriteResults2.open("results_2.dat");
    if (WriteResults2.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults2 << sumP_prog[i] << " " << errP_prog[i] << " " << "\t" << endl;
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
