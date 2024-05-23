//
// Created by francesco on 23/05/24.
//

#include "Population.h"

using namespace std;
using namespace arma;

Population::Population() {}
// Default constructor, does not perform any action

Population::~Population() {}
// Default destructor, does not perform any action

Path Population::getPath(int i) {
    return population(i);
}
