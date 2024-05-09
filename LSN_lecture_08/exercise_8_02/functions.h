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

double potential(double x);

double kinetic_energy(double x, double mu, double sigma);

double psi_trial(double x, double mu, double sigma);

double pdf_function(double, double mu, double sigma);

//Metropolis algorithm with Uniform extraction
std::tuple<std::vector<double>, std::vector<double>>
Metropolis_Uniform(double x, Random rnd, double metropolis_step,
                   double (*pdf_function)(double x, double mu, double sigma),
                   double (*potential)(double x),
                   double (*kinetic_energy)(double x, double mu, double sigma), double mu, double sigma, int blocks,
                   int steps,
                   std::string filename);

//Cumulative average and error
void cumulativeAverage(std::vector<double> average, std::vector<double> average2, std::string filename);

//Error
double error(std::vector<double> av, std::vector<double> av2, int n);



#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
