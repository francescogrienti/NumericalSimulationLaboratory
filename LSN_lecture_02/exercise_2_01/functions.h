//
// Created by francesco on 17/03/24.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

#ifndef __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
#define __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__

//Mean
std::tuple<std::vector<double>, std::vector<double>> mean(int M, int N, Random rnd);

//Mean of the standard deviation
std::tuple<std::vector<double>, std::vector<double>> mean(int M, int N, Random rnd, double mu);

//Cumulative average and error
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
cumulativeAverage(std::vector<double> average, std::vector<double> average2);

//Writing on file
void writeOnFile(std::vector<double> average, std::vector<double> error, std::string filename);

void writeOnFile(std::vector<double> chi2, std::string filename);

//Exponential writing to file
void writeOnFile(int M, std::vector<int> N, Random rnd, double lambda,
                 std::vector<double> (*function)(Random, int, int, double), std::string filename);

//Uniform writing to file
void writeOnFile(int M, std::vector<int> N, Random rnd, std::vector<double> (*function)(Random, int, int),
                 std::string filename);

//Cauchy writing to file
void writeOnFile(int M, std::vector<int> N, Random rnd, double mean, double gamma,
                 std::vector<double> (*function)(Random, int, int, double, double),
                 std::string filename);

//Error
double error(std::vector<double> av, std::vector<double> av2, int n);

//Integrand function
double func_f(double x);

//Probability function
double prob_p(double x);

//Auxiliary function
double func_g(double x);

//Initialize random number generator
Random initialize(Random rnd, std::vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file);


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
