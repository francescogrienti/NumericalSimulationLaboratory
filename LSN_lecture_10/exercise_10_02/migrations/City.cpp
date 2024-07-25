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

double City::getCoordinate() {
    return angular_coordinate;
}

std::vector<double> City::getSquareCoordinates() {
    vector<double> coordinates = {x, y};
    return coordinates;
}

//Setter methods
void City::setLabel(int label) {
    city_label = label;
    return;
}

void City::setCoordinate(double coordinate) {
    angular_coordinate = coordinate;
    return;
}

void City::setSquareCoordinates(double x_1, double x_2) {
    this->x = x_1;
    this->y = x_2;
    return;
}

