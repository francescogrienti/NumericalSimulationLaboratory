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
        for (int i = 0; i < (int)(average.size()); i++) {
            WriteResults << average[i] << " " << error[i] << " " << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

void writeOnFile(std::vector<double> chi2, std::string filename){
    ofstream WriteResults;
    WriteResults.open(filename);
    if (WriteResults.is_open()) {
        for (int i = 0; i < (int)(chi2.size()); i++) {
            WriteResults << chi2[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults.close();
}

tuple<vector<double>, vector<double>, vector<double>>
cumulativeAverage(vector<double> average, vector<double> average2) {
    vector<double> sumProg((int)(average.size()), 0.);
    vector<double> sum2Prog((int)(average.size()), 0.);
    vector<double> errProg((int)(average.size()), 0.);
    for (int k = 0; k < (int)(average.size()); k++) {
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


double error(vector<double> av, vector<double> av2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n], 2)) / (sqrt(n));
    }
}




