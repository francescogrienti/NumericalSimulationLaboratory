//
// Created by francesco on 07/06/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_GENETICS_H
#define NUMERICALSIMULATIONLABORATORY_GENETICS_H

#include "random.h"
#include <unordered_map>
#include <unordered_set>
#include "Path.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>   // per rand() e srand()
#include <ctime>     // per time()
#include <utility>   // per std::pair

class Genetics {
private:
    Path path;
    int n_cities;
    int pop_size;
public:
    // Default constructor
    Genetics();

    // Destructor
    ~Genetics();

    //Set number of cities in the path
    void setCitiesPath(int n);

    //Set number of individuals in the population
    void setPopSize(int n);

    //Initialize the map
    void initialize_path(Random &rnd, string type);

    //Creation of the first population
    vector<vector<int>> first_pop(Random &rnd);

    //Function checking that the individual fulfils the bonds
    bool check_function(std::vector<int> &labels);

    //Function sorting the population
    void sort_paths(std::vector<vector<int>> &population, string type);

    //Selection function
    int selection_operator(vector<vector<int>> &population, Random &rnd, int p);

    //Pair permutation operator
    void pair_permutation(double prob, vector<int> &labels, Random &rnd);

    //Shift operator
    void shift_operator(double prob, vector<int> &labels, int N_elem, int shift, Random &rnd);

    //M-permutation
    void m_permutation(double prob, vector<int> &labels, Random &rnd);

    //Inverse operator
    void inverse_operator(double prob, vector<int> &labels, Random &rnd);

    //Computing distances
    double compute_best_path(vector<int> &labels, double r, string type);

    //Computing average of the best half of the population
    double compute_half_best_path(vector<vector<int>> &population, double r, string type);

    //Cross-over operator
    pair<vector<int>, vector<int>> cross_over_operator(vector<int> &parent_1, vector<int> &parent_2, Random &rnd);

    //Get coordinates of the city
    double getCityCoordinate(int i);

    //Get coordinates of the city
    std::vector<double> getCitySquareCoordinates(int i);
};


#endif //NUMERICALSIMULATIONLABORATORY_GENETICS_H
