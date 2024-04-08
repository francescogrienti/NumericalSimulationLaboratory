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


#endif //__NUMERICALSIMULATIONLABORATORY_FUNCTIONS_H__
