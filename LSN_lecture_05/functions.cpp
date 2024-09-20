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
            WriteResults << average[i] << " " << error[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}


void cumulativeAverage(vector<double> average, vector<double> average2, string filename) {
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

    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < (int) (average.size()); i++) {
            WriteResults << sumProg[i] << " " << errProg[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
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


double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}

double acceptance(vector<double> x, vector<double> x_n) {
    //Add probability density function p(x_n)/p(x)
    return min(1., x_n[0] / x[0]);
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

double pdf_wave_function_GS(std::vector<double> x) {
    double a_0 = 1.;
    return (pow(a_0, (-3. / 2.)) / sqrt(M_PI)) * exp((-1) * (sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) / a_0);
}

double pdf_wave_function_ES(std::vector<double> x) {
    double a_0 = 1.;
    return (pow(a_0, (-5. / 2.)) / (8. * sqrt(M_PI))) * sqrt(2.) * x[2] *
           exp((-1) * (sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))) / (2. * a_0));
}

std::tuple<vector<double>, vector<double>>
Metropolis_Gauss(vector<double> x, Random rnd, double metropolis_step, double (*pdf_function)(vector<double>),
                 int blocks, int steps, string filename) {

    double radius_sum = 0.;
    double acceptance = 0.;
    vector<double> ave(blocks, 0.);
    vector<double> ave2(blocks, 0.);
    double r = 0.;
    vector<double> x_k(3, 0.);
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < blocks; i++) {
            radius_sum = 0.;
            for (int j = 0; j < (steps / blocks); j++) {
                x_k[0] = rnd.Gauss(x[0], metropolis_step);
                x_k[1] = rnd.Gauss(x[1], metropolis_step);
                x_k[2] = rnd.Gauss(x[2], metropolis_step);
                acceptance = min(1., pow(pdf_function(x_k), 2) / pow(pdf_function(x), 2));
                r = rnd.Rannyu();
                if (r <= acceptance) {
                    x[0] = x_k[0];
                    x[1] = x_k[1];
                    x[2] = x_k[2];
                }
                radius_sum += sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
                WriteResults << acceptance << " " << x[0] << " " << x[1] << " " << x[2] << " " << "\t" << endl;
            }
            ave[i] = radius_sum / (steps / blocks);
            ave2[i] = double(pow(ave[i], 2));
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
    return make_tuple(ave, ave2);

}

std::tuple<vector<double>, vector<double>>
Metropolis_Uniform(vector<double> x, Random rnd, double metropolis_step, double (*pdf_function)(vector<double>),
                   int blocks, int steps, string filename) {

    double radius_sum = 0.;
    double acceptance = 0.;
    vector<double> ave(blocks, 0.);
    vector<double> ave2(blocks, 0.);
    double r = 0.;
    vector<double> x_k(3, 0.);
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < blocks; i++) {
            radius_sum = 0.;
            for (int j = 0; j < (steps / blocks); j++) {
                x_k[0] = rnd.Rannyu(x[0] - metropolis_step, x[0] + metropolis_step);
                x_k[1] = rnd.Rannyu(x[1] - metropolis_step, x[1] + metropolis_step);
                x_k[2] = rnd.Rannyu(x[2] - metropolis_step, x[2] + metropolis_step);
                acceptance = min(1., pow(pdf_function(x_k), 2) / pow(pdf_function(x), 2));
                r = rnd.Rannyu();
                if (r <= acceptance) {
                    x[0] = x_k[0];
                    x[1] = x_k[1];
                    x[2] = x_k[2];
                }
                radius_sum += sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
                WriteResults << acceptance << " " << x[0] << " " << x[1] << " " << x[2] << " " << "\t" << endl;
            }
            ave[i] = radius_sum / (steps / blocks);
            ave2[i] = double(pow(ave[i], 2));
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
    return make_tuple(ave, ave2);

}

void Metropolis_Uniform_eq(vector<double> x, Random rnd, double metropolis_step, double (*pdf_function)(vector<double>),
                   int blocks, int steps, string filename) {

    double radius = 0.;
    double acceptance = 0.;
    double r = 0.;
    vector<double> x_k(3, 0.);
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < blocks; i++) {
            radius = 0.;
            for (int j = 0; j < (steps / blocks); j++) {
                x_k[0] = rnd.Rannyu(x[0] - metropolis_step, x[0] + metropolis_step);
                x_k[1] = rnd.Rannyu(x[1] - metropolis_step, x[1] + metropolis_step);
                x_k[2] = rnd.Rannyu(x[2] - metropolis_step, x[2] + metropolis_step);
                acceptance = min(1., pow(pdf_function(x_k), 2) / pow(pdf_function(x), 2));
                r = rnd.Rannyu();
                if (r <= acceptance) {
                    x[0] = x_k[0];
                    x[1] = x_k[1];
                    x[2] = x_k[2];
                }
                radius = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
                WriteResults << radius << " " << x[0] << " " << x[1] << " " << x[2] << " " << "\t" << endl;
            }
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
    return;

}

