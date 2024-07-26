//
// Created by francesco on 23/05/24.
//

#include <cstdlib>
#include "City.h"
#include <vector>


using namespace std;

City::City() {}
// Default constructor, does not perform any action

City::~City() {}
// Default destructor, does not perform any action

//Getter methods
int City::getLabel() {
    return city_label;
}


std::vector<double> City::getProvinceCoordinates() {
    vector<double> coordinates = {x, y};
    return coordinates;
}

//Setter methods
void City::setLabel(int label) {
    city_label = label;
    return;
}

void City::setProvinceCoordinates(double x_1, double x_2) {
    this->x = x_1;
    this->y = x_2;
    return;
}

