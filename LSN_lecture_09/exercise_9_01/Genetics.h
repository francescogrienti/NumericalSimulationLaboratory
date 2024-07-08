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
    vector<double> probabilities;
public:
    // Default constructor
    Genetics();

    // Destructor
    ~Genetics();

    //Set the probability for the mutation operators
    void setProbabilities(std::vector<double> prob);

    //Get the probability for the mutation operators
    vector<double> getProbabilities();

    //Set number of cities in the path
    void setCitiesPath(int n);

    //Set number of individuals in the population
    void setPopSize(int n);

    //Initialize the map -- circle
    void initialize_path_circle(Random &rnd);

    //Initialize the map -- square
    void initialize_path_square(Random &rnd);

    //Creation of the first population
    vector<vector<int>> first_pop(Random &rnd);

    //Function checking that the individual fulfils the bonds
    bool check_function(std::vector<int> &labels);

    //Function sorting the population
    void sort_paths(std::vector<vector<int>> &population);

    //Selection function
    vector<int> selection_operator(const vector<vector<int>> &population, Random &rnd, int p);

    //Pair permutation operator
    void pair_permutation(double prob, vector<int> &labels, Random &rnd);

    //Shift operator
    void shift_operator(double prob, vector<int> &labels, int N_elem, int shift, Random &rnd);

    //M-permutation
    void m_permutation(double prob, vector<int> &labels, int n, Random &rnd);

    //Inverse operator
    void inverse_operator(double prob, vector<int> &labels, int n, Random &rnd);

    //Function for the mutation
    void mutation(Random &rnd);

    //Computing distances
    double compute_best_path(vector<int> &labels, double r);

    //Computing average of the best half of the population
    double compute_half_best_path(vector<vector<int>> &population, double r);

    //Cross-over operator
    pair<vector<int>, vector<int>> cross_over_operator(vector<int> &parent_1, vector<int> &parent_2, Random &rnd);

    //Get coordinates of the city
    double getCityCoordinate(int i);

};


#endif //NUMERICALSIMULATIONLABORATORY_GENETICS_H
