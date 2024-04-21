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

//Initialize random number generator
Random initialize(Random rnd, std::vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
sampling(Random rnd, double r, double sigma, int blocks, int steps, int S_0, int T, int strike_price, double t_i,
         std::string sampling_type);

//Cumulative average and error
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
cumulativeAverage(std::vector<double> average, std::vector<double> average2);

//Error
double error(std::vector<double> av, std::vector<double> av2, int n);

//Writing on file
void writeOnFile(std::vector<double> average, std::vector<double> error, std::string filename);


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
