
//
// Created by francesco on 17/03/24.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include "functions.h"
#include "random.h"

using namespace std;

Random initialize(Random rnd, vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file) {
    ifstream Primes(prime_file);
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input(input_file);
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

    return rnd;
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>>
sampling(Random rnd, double r, double sigma, int blocks, int steps, int S_0, int T, int strike_price, double t_i,
         string sampling_type) {

    vector<double> ave_C_i(blocks, 0.);
    vector<double> ave2_C_i(blocks, 0.);
    vector<double> ave_P_i(blocks, 0.);
    vector<double> ave2_P_i(blocks, 0.);

    if (sampling_type == "discretized") {
        vector<double> S(blocks + 1, 0.);
        for (int i = 0; i < blocks; i++) {
            double sum_C = 0.;
            double sum_P = 0.;
            for (int k = 0; k < (steps / blocks); k++) {
                double z_i = rnd.Gauss(0., 1.);
                S[0] = S_0;
                S[1] = S_0 * exp(((r - pow(sigma, 2) / 2)) * t_i + sigma * z_i * sqrt(t_i));
                for (int j = 1; j < blocks; j++) {
                    z_i = rnd.Gauss(0., 1.);
                    S[j + 1] = S[j] * exp(((r - pow(sigma, 2) / 2)) * t_i + sigma * z_i * sqrt(t_i));
                }
                sum_C += exp((-1) * r * T) * max(0., S[blocks] - strike_price);
                sum_P += exp((-1) * r * T) * max(0., strike_price - S[blocks]);
            }
            ave_C_i[i] = sum_C / (steps / blocks);
            ave_P_i[i] = sum_P / (steps / blocks);
            ave2_C_i[i] = pow(ave_C_i[i], 2);
            ave2_P_i[i] = pow(ave_P_i[i], 2);
        }
    } else if (sampling_type == "direct") {
        double S_T = 0.;
        for (int i = 0; i < blocks; i++) {
            double sum_C = 0.;
            double sum_P = 0.;
            for (int k = 0; k < (steps / blocks); k++) {
                double z_i = rnd.Gauss(0., 1.);
                S_T = S_0 * exp(((r - pow(sigma, 2) / 2)) * T + sigma * z_i * sqrt(T));
                sum_C += exp((-1) * r * T) * max(0., S_T - strike_price);
                sum_P += exp((-1) * r * T) * max(0., strike_price - S_T);
            }
            ave_C_i[i] = sum_C / (steps / blocks);
            ave_P_i[i] = sum_P / (steps / blocks);
            ave2_C_i[i] = pow(ave_C_i[i], 2);
            ave2_P_i[i] = pow(ave_P_i[i], 2);
        }
    }
    return make_tuple(ave_C_i, ave2_C_i, ave_P_i, ave2_P_i);
}

tuple<vector<double>, vector<double>, vector<double>>
cumulativeAverage(vector<double> average, vector<double> average2) {
    vector<double> sumProg((int) (average.size()), 0.);
    vector<double> sum2Prog((int) (average.size()), 0.);
    vector<double> errProg((int) (average.size()), 0.);
    for (int k = 0; k < (int) (average.size()); k++) {
        sumProg[k] = 0.;
        sum2Prog[k] = 0.;
        for (int l = 0; l < k + 1; l++) {
            sumProg[k] += average[l];
            sum2Prog[k] += average2[l];
        }
        sumProg[k] /= (k + 1); //Cumulative average of the standard deviation
        sum2Prog[k] /= (k + 1); //Cumulative square average of the standard deviation
        errProg[k] = error(sumProg, sum2Prog, k); //Statistical uncertainty
    }

    return make_tuple(sumProg, sum2Prog, errProg);
}

double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}

void writeOnFile(vector<double> average, vector<double> error, string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < (int) (average.size()); i++) {
            WriteResults << average[i] << " " << error[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}







