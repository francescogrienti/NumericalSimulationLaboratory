//
// Created by francesco on 23/05/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_PATH_H
#define NUMERICALSIMULATIONLABORATORY_PATH_H

#include <armadillo>
#include "City.h"
#include "functions.h"
#include <vector>

using namespace std;
using namespace arma;

class Path {
private:
    field<City> path;
    double path_length;
public:
    // Default constructor
    Path();

    // Destructor
    ~Path();

    void initialize(const int n_cities);

    //Method for getting the city in the path
    City getCity(int i);

    //Method for getting the length of the path
    double getPathLength();

    //Method for setting the length
    void setPathLength(double length);

    //Method for setting the vector of cities in the path
    void setCity(int label, double coordinate, int i);

    //Method for getting the vector of labels
    std::vector<int> getLabels();

};


#endif //NUMERICALSIMULATIONLABORATORY_PATH_H
