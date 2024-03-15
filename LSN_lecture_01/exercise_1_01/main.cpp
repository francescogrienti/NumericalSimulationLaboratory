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
#include "random.h"

using namespace std;

double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}

int main(int argc, char *argv[]) {

    Random rnd;
    int M = 100000; //Number of total throws
    int N = 100; //Number of blocks
    int L = M / N; //Number of throws per block
    double mu = 0.5;
    vector<double> ave(N, 0.); //Arrays for storing the average and the square of the average for each block
    vector<double> ave2(N, 0.); //Arrays for storing the average and the square of the average for each block
    vector<double> sumProg(N, 0.);
    vector<double> sum2Prog(N, 0.);
    vector<double> errProg(N, 0.);
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
        double sum = 0.;
        for (int j = 0; j < L; j++) {
            double r = rnd.Rannyu();
            sum += r;
        }
        ave[i] = double(sum / L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
    }
    for (int k = 0; k < N; k++) {
        sumProg[k] = 0.;
        sum2Prog[k] = 0.;
        for (int l = 0; l < k + 1; l++) {
            sumProg[k] += ave[l];
            sum2Prog[k] += ave2[l];
        }
        sumProg[k] /= (k + 1); //Cumulative average
        sum2Prog[k] /= (k + 1); //Cumulative square average
        errProg[k] = error(sumProg, sum2Prog, k); //Statistical uncertainty
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults1 << sumProg[i] << " " << errProg[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults1.close();

    //Block of code for the evaluation of the standard deviation
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        for (int j = 0; j < L; j++) {
            double r = rnd.Rannyu();
            sum += pow(r - mu, 2);
        }
        ave[i] = sum / L; //Store average values of the standard deviation for each block
        ave2[i] = pow(ave[i], 2); //Store square of the average of the standard deviation for each block
    }

    for (int k = 0; k < N; k++) {
        sumProg[k] = 0.;
        sum2Prog[k] = 0.;
        for (int l = 0; l < k + 1; l++) {
            sumProg[k] += ave[l];
            sum2Prog[k] += ave2[l];
        }
        sumProg[k] /= (k + 1); //Cumulative average of the standard deviation
        sum2Prog[k] /= (k + 1); //Cumulative square average of the standard deviation
        errProg[k] = error(sumProg, sum2Prog, k); //Statistical uncertainty
    }

    ofstream WriteResults2;
    WriteResults2.open("results_2.dat");
    if (WriteResults2.is_open()) {
        for (int i = 0; i < N; i++) {
            WriteResults2 << sumProg[i] << " " << errProg[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults2.close();

    //Chi-2
    int subI = 100; //Number of sub-intervals of [0,1)
    int n = 10000; //Number of throws
    int expecValue = n / subI;
    vector<double> chi2(subI, 0.);
    vector<int> count(subI, 0);

    for (int i = 0; i < subI; i++) {
        for (int k = 0; k < n; k++) {
            double r = rnd.Rannyu();
            count[(int) (r * subI)] += 1;
        }
        for (int j = 0; j < subI; j++) {
            chi2[j] += double(pow(count[j] - expecValue, 2)) / (double) (expecValue);
            count[j] = 0;
        }
    }

    ofstream WriteResults3;
    WriteResults3.open("results_3.dat");
    if (WriteResults3.is_open()) {
        for (int i = 0; i < subI; i++) {
            WriteResults3 << chi2[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults3.close();

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
