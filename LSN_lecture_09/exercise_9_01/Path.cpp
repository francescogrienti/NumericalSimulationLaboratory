//
// Created by francesco on 23/05/24.
//

#include "Path.h"
#include "City.h"
#include <cmath>
#include <string>

using namespace std;

Path::Path() {}
// Default constructor, does not perform any action

Path::~Path() {}
// Default destructor, does not perform any action

//Initialize path

void Path::initialize_path(int n) {
    path.resize(n);
    return;
}

//Getter methods
City Path::getCity(int i) {
    return path[i - 1];
}


void Path::setCity(int label, double coordinate, int i) {
    path[i].setLabel(label);
    path[i].setCoordinate(coordinate);
    return;
}

void Path::setCity(int label, double x_1, double x_2, int i) {
    path[i].setLabel(label);
    path[i].setSquareCoordinates(x_1, x_2);
    return;
}

std::vector<int> Path::getLabels() {
    int size = path.size();
    vector<int> labels(size, 0);
    for (int i = 0; i < size; i++) {
        labels[i] = path[i].getLabel();
    }
    return labels;
}

//Function for computing the distance between cities (using L-1 norm)
double Path::L1_norm(City city_1, City city_2, double r, string type) {
    vector<double> x(2, 0.);
    vector<double> y(2, 0.);
    if (type == "circle") {
        x[0] = r * cos(city_1.getCoordinate());
        x[1] = r * sin(city_1.getCoordinate());
        y[0] = r * cos(city_2.getCoordinate());
        y[1] = r * sin(city_2.getCoordinate());
    } else if (type == "square") {
        x[0] = city_1.getSquareCoordinates()[0]; //X City 1
        x[1] = city_1.getSquareCoordinates()[1]; //Y City 1
        y[0] = city_2.getSquareCoordinates()[0]; //X City 2
        y[1] = city_2.getSquareCoordinates()[1]; //Y City 2
    }
    return sqrt((pow(x[0] - y[0], 2)) + (pow(x[1] - y[1], 2)));
}

