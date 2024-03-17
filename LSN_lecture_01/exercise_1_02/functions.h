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

//Function returning the sum of n random variables sampled from an exponential distribution
std::vector<double> expon_prob_sum(Random rnd, int n, int throws, double lambda);

//Function returning the sum of n random variables sampled from a uniform distribution
std::vector<double> uniform_prob_sum(Random rnd, int n, int throws);

//Function returning the sum of n random variables sampled from a uniform distribution
std::vector<double> cauchy_prob_sum(Random rnd, int n, int throws, double mean, double gamma);


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
