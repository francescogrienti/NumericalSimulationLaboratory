
//
// Created by francesco on 17/03/24.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <cmath>
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

double L1_norm(City city1, City city2, double r) {
    vector<double> x(2, 0.);
    vector<double> y(2, 0.);
    x[0] = r * cos(city1.getCoordinate());
    x[1] = r * sin(city1.getCoordinate());
    y[0] = r * cos(city2.getCoordinate());
    y[1] = r * sin(city2.getCoordinate());
    return sqrt((pow(x[0] - y[0], 2)) + (pow(x[1] - y[1], 2)));
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

//Check function
bool check_function(std::vector<int> &labels) {
    std::unordered_map<int, bool> dict;
    for (int i = 1; i < labels.size() - 1; i++) {
        if (dict.find(i) != dict.end()) {
            return false;
        }
        dict[i] = true;
    }
    return true;
}

//Selection operator
//Path selection_operator(field<Path> population)

double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}









