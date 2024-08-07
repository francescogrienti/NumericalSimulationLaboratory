//
// Created by francesco on 23/05/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_PATH_H
#define NUMERICALSIMULATIONLABORATORY_PATH_H

#include "City.h"
#include <vector>
#include <string>

using namespace std;

class Path {
private:
    vector<City> path;
public:
    // Default constructor
    Path();

    // Destructor
    ~Path();

    void initialize_path(int n);

    //Method for getting the city in the path
    City getCity(int i);

    //Function for computing the distance between cities (using L-1 norm)
    double L1_norm(City city_1, City city_2);

    //Method for setting the vector of cities in the path (Square)
    void setCity(int label, double x_1, double x_2, int i);

    //Method for getting the vector of labels
    std::vector<int> getLabels();
};


#endif //NUMERICALSIMULATIONLABORATORY_PATH_H
