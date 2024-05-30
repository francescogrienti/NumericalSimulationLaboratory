//
// Created by francesco on 17/03/24.
//
#include <vector>
#include <unordered_map>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "City.h"
#include "Path.h"

#ifndef __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
#define __NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__

//Initialize random number generator
Random initialize(Random rnd, std::vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file);

//Function for computing the distance between cities (using L-1 norm)
double L1_norm(City city_1, City city_2, double r);

//Function checking that the individual fulfils the bonds
bool check_function(std::vector<int>& labels);

//Cumulative average and error
void cumulativeAverage(std::vector<double> average, std::vector<double> average2, std::string filename);

//Error
double error(std::vector<double> av, std::vector<double> av2, int n);


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
