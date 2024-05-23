//
// Created by francesco on 23/05/24.
//

#include "Path.h"
#include "City.h"

using namespace std;
using namespace arma;

Path::Path() {}
// Default constructor, does not perform any action

Path::~Path() {}
// Default destructor, does not perform any action

//Getter methods
City Path::getCity(int i) {
    return path(i);
}

double Path::getPathLength() {
    return path_length;
}

//Setter methods
void Path::setPathLength(double length) {
    path_length = length;
    return;
}

void Path::setCity(int label, double coordinate, int i) {
    path(i).setLabel(label);
    path(i).setCoordinate(coordinate);
    return;
}