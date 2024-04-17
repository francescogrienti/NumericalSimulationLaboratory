//
// Created by francesco on 17/03/24.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "algorithm"


#ifndef __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
#define __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__

//Mean
std::tuple<std::vector<double>, std::vector<double>> mean(int M, int N, Random rnd);

//Mean of the standard deviation
std::tuple<std::vector<double>, std::vector<double>> mean(int M, int N, Random rnd, double mu);

//Cumulative average and error
void cumulativeAverage(std::vector<double> average, std::vector<double> average2, std::string filename);

//Writing on file
void writeOnFile(std::vector<double> average, std::vector<double> error, std::string filename);

//Error
double error(std::vector<double> av, std::vector<double> av2, int n);

//Acceptance function
double acceptance(std::vector<double> x, std::vector<double> x_n);

//Metropolis algorithm with Gaussian extraction
std::tuple<std::vector<double>, std::vector<double>>
Metropolis_Gauss(std::vector<double> x, Random rnd, double metropolis_step,
                 double (*pdf_function)(std::vector<double>), int blocks, int steps, std::string filename);

//Metropolis algorithm with Uniform extraction
std::tuple<std::vector<double>, std::vector<double>>
Metropolis_Uniform(std::vector<double> x, Random rnd, double metropolis_step,
                   double (*pdf_function)(std::vector<double>), int blocks, int steps, std::string filename);

//Probability density function for the ground state
double pdf_wave_function_GS(std::vector<double> x);

//Probability density function for the excited state
double pdf_wave_function_ES(std::vector<double> x);

//

//Initialize random number generator
Random initialize(Random rnd, std::vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file);


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
