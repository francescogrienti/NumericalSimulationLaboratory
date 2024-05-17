
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

double potential(double x) {
    return pow(x, 4) - (5. / 2.) * (pow(x, 2));
}

double psi_trial(double x, double mu, double sigma) {
    return exp((-1.) * (pow(x - mu, 2) / (2 * pow(sigma, 2)))) + exp((-1.) * (pow(x + mu, 2) / (2 * pow(sigma, 2))));
}

double kinetic_energy(double x, double mu, double sigma) {
    double first_term = exp(-(pow(x - mu, 2) / (2 * pow(sigma, 2)))) *
                        ((pow(x, 2) - 2 * x * mu + pow(mu, 2)) / pow(sigma, 4) - 1 / pow(sigma, 2));

    double second_term = exp(-(pow(x + mu, 2) / (2 * pow(sigma, 2)))) *
                         ((pow(x, 2) + 2 * x * mu + pow(mu, 2)) / pow(sigma, 4) - 1 / pow(sigma, 2));

    return (-0.5) * ((first_term + second_term) /
                     (exp((-1.) * (pow(x - mu, 2) / (2 * pow(sigma, 2)))) +
                      exp((-1.) * (pow(x + mu, 2) / (2 * pow(sigma, 2))))));

}

std::tuple<std::vector<double>, std::vector<double>>
Metropolis_Uniform(double x, Random rnd, double metropolis_step,
                   double (*pdf_function)(double x, double mu, double sigma),
                   double (*potential)(double x),
                   double (*kinetic_energy)(double x, double mu, double sigma), double mu, double sigma, int blocks,
                   int steps,
                   std::string filename) {

    double integral_sum = 0.;
    double acceptance = 0.;
    vector<double> ave(blocks, 0.);
    vector<double> ave2(blocks, 0.);
    double r = 0.;
    double x_k = 0.;
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < blocks; i++) {
            integral_sum = 0.;
            for (int j = 0; j < (steps / blocks); j++) {
                x_k = rnd.Rannyu(x - metropolis_step, x + metropolis_step);
                acceptance = min(1., pow(pdf_function(x_k, mu, sigma), 2) /
                                     pow(pdf_function(x, mu, sigma), 2));
                r = rnd.Rannyu();
                if (r <= acceptance) {
                    x = x_k;
                }
                integral_sum += kinetic_energy(x, mu, sigma) + potential(x);
                WriteResults << acceptance << " " << x << " " << "\t" << endl;
            }
            ave[i] = integral_sum / (steps / blocks);
            ave2[i] = double(pow(ave[i], 2));
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
    return make_tuple(ave, ave2);
}

double pdf_function(double x, double mu, double sigma) {
    double psi_trial_2 = pow(
            exp((-1.) * (pow(x - mu, 2) / (2 * pow(sigma, 2)))) + exp((-1.) * (pow(x + mu, 2) / (2 * pow(sigma, 2)))),
            2);
    return psi_trial_2;
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

double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}







