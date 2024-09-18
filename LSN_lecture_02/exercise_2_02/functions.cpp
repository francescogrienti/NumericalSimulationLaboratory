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

void writeOnFile(vector<double> average, vector<double> error, string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < (int) (average.size()); i++) {
            WriteResults << sqrt(average[i]) << " " << error[i] / (2 * sqrt(average[i])) << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

void writeOnFile(std::vector<double> chi2, std::string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < (int) (chi2.size()); i++) {
            WriteResults << chi2[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

void writeOnFile(int M, std::vector<int> N, Random rnd, double lambda,
                 std::vector<double> (*function)(Random, int, int, double), std::string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    vector<double> sum1 = function(rnd, N[0], M, lambda);
    vector<double> sum2 = function(rnd, N[1], M, lambda);
    vector<double> sum10 = function(rnd, N[2], M, lambda);
    vector<double> sum100 = function(rnd, N[3], M, lambda);
    if (WriteResults.is_open()) {
        for (int i = 0; i < M; i++) {
            WriteResults << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

void writeOnFile(int M, std::vector<int> N, Random rnd,
                 std::vector<double> (*function)(Random, int, int), std::string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    vector<double> sum1 = function(rnd, N[0], M);
    vector<double> sum2 = function(rnd, N[1], M);
    vector<double> sum10 = function(rnd, N[2], M);
    vector<double> sum100 = function(rnd, N[3], M);
    if (WriteResults.is_open()) {
        for (int i = 0; i < M; i++) {
            WriteResults << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

void writeOnFile(int M, std::vector<int> N, Random rnd, double mean, double gamma,
                 std::vector<double> (*function)(Random, int, int, double, double), std::string filename) {
    ofstream WriteResults;
    WriteResults.open(filename);
    vector<double> sum1 = function(rnd, N[0], M, mean, gamma);
    vector<double> sum2 = function(rnd, N[1], M, mean, gamma);
    vector<double> sum10 = function(rnd, N[2], M, mean, gamma);
    vector<double> sum100 = function(rnd, N[3], M, mean, gamma);
    if (WriteResults.is_open()) {
        for (int i = 0; i < M; i++) {
            WriteResults << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
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

tuple<vector<double>, vector<double>> mean(int M, int N, Random rnd) {
    int L = M / N;
    vector<double> ave(N, 0.);
    vector<double> ave2(N, 0.);
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        for (int j = 0; j < L; j++) {
            double r = rnd.Rannyu();
            sum += r;
        }
        ave[i] = double(sum / L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square-average for each block
    }

    return make_tuple(ave, ave2);
}

tuple<vector<double>, vector<double>> mean(int M, int N, Random rnd, double mu) {
    int L = M / N;
    vector<double> ave(N, 0.);
    vector<double> ave2(N, 0.);
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        for (int j = 0; j < L; j++) {
            double r = rnd.Rannyu();
            sum += pow(r - mu, 2);
        }
        ave[i] = double(sum / L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square-average for each block
    }

    return make_tuple(ave, ave2);
}

tuple<vector<double>, vector<double>> mean_and_error(vector<vector<double>> matrix_case, int N, int L) {
    vector<vector<double>> ave(N, vector<double>(L, 0.));
    vector<vector<double>> ave2(N, vector<double>(L, 0.));
    vector<double> average(N, 0.);
    vector<double> error(N, 0.);
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < L; i++) {
            ave[k][i] = (matrix_case[k][i] / L);
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

    return make_tuple(average, error);
}


double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}

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


double func_f(double x) {
    return (M_PI / 2.) * cos((M_PI / 2.) * x);
}

double prob_p(double x) {
    return (1.5) * (1 - pow(x, 2));
}

double func_g(double x) {
    return ((M_PI / 3.) * cos((M_PI / 2.) * x)) / (1 - pow(x, 2));
}






