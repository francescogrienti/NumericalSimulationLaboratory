//
// Created by francesco on 23/05/24.
//

#include "Path.h"
#include "City.h"
#include <cmath>

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


void Path::setCity(int label, double x_1, double x_2, int i) {
    path[i].setLabel(label);
    path[i].setProvinceCoordinates(x_1, x_2);
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
double Path::L1_norm(City city_1, City city_2) {
    vector<double> x(2, 0.);
    vector<double> y(2, 0.);
    x[0] = city_1.getProvinceCoordinates()[0]; //X City 1
    x[1] = city_1.getProvinceCoordinates()[1]; //Y City 1
    y[0] = city_2.getProvinceCoordinates()[0]; //X City 2
    y[1] = city_2.getProvinceCoordinates()[1]; //Y City 2

    return sqrt((pow(x[0] - y[0], 2)) + (pow(x[1] - y[1], 2)));
}

