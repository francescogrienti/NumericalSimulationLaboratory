//
// Created by francesco on 23/05/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_POPULATION_H
#define NUMERICALSIMULATIONLABORATORY_POPULATION_H

#include "Path.h"
#include <vector>

class Population {
private:
    field<Path> population;
public:
    // Default constructor
    Population();

    // Destructor
    ~Population();

    //Method for getting i-th path in the population
    Path getPath(int i);
};


#endif //NUMERICALSIMULATIONLABORATORY_POPULATION_H
