//
// Created by francesco on 07/06/24.
//

#include "Genetics.h"
#include <vector>
#include <tuple>
#include <cmath>
#include "random.h"
#include <algorithm>
#include <string>

using namespace std;

Genetics::Genetics() {}
// Default constructor, does not perform any action

Genetics::~Genetics() {}
// Default destructor, does not perform any action

void Genetics::setCitiesPath(int n) {
    n_cities = n;
    return;
}

void Genetics::setPopSize(int n) {
    pop_size = n;
    return;
}

//Initialize the map
void Genetics::initialize_path(string filename) {
    vector<vector<double>> coordinates;
    vector<int> labels(n_cities, 0);
    for (int k = 0; k < n_cities; k++) {
        labels[k] = k + 1;
    }

    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Unable to open file";
        return; // return with error code
    }

    double value1, value2;
    while (inputFile >> value1 >> value2) {
        std::vector<double> row = {value1, value2};
        coordinates.push_back(row);
    }
    inputFile.close();
    path.initialize_path(n_cities);
    for (int k = 0; k < n_cities; k++) {
        path.setCity(labels[k], coordinates[k][0], coordinates[k][1], k);
    }
    return;
}

//Creation of the first population using the mutation operator
vector<vector<int>> Genetics::first_pop(Random &rnd) {
    vector<vector<int>> population(pop_size, vector<int>(n_cities + 1, 0));

    //GENERATION OF THE POPULATION
    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j <= n_cities; j++) {
            population[i][0] = 1;
            population[i][n_cities] = 1;
            population[i][j + 1] = j + 2;
        }
    }
    for (int i = 0; i < pop_size; i++) {
        random_shuffle(population[i].begin() + 1, population[i].end() - 1);
    }
    return population;
}

bool Genetics::check_function(std::vector<int> &labels) {
    std::unordered_map<int, bool> dict;
    for (int i = 1; i < labels.size() - 1; i++) {
        if (dict.find(i) != dict.end()) {
            return false;
        }
        dict[i] = true;
    }
    return true;
}

void Genetics::sort_paths(std::vector<vector<int>> &population) {
    vector<double> path_length(population.size(), 0.0);

    // Compute the distances of the various paths
    for (int i = 0; i < population.size(); i++) {
        for (int j = 0; j < n_cities; j++) {
            path_length[i] += path.L1_norm(path.getCity(population[i][j]), path.getCity(population[i][j + 1]));
        }
    }

    // Sorting
    // Create the struct to use for sorting
    struct VecWithVal {
        std::vector<int> vec;
        double value;
    };

    std::vector<VecWithVal> vecWithVals(population.size());
    for (int i = 0; i < population.size(); i++) {
        vecWithVals[i].vec = population[i];
        vecWithVals[i].value = path_length[i];
    }

    auto compare = [](const VecWithVal &a, const VecWithVal &b) {
        return a.value < b.value;
    };

    std::sort(vecWithVals.begin(), vecWithVals.end(), compare);

    for (int i = 0; i < population.size(); i++) {
        population[i] = vecWithVals[i].vec;
    }

    return;
}

void Genetics::pair_permutation(double prob, vector<int> &labels, Random &rnd) {
    if (rnd.Rannyu() < prob) {
        //Generation of two random integer numbers to select the pair indeces to be muted
        int n_1 = int(rnd.Rannyu(1., n_cities));
        int n_2 = int(rnd.Rannyu(1., n_cities));
        while (n_1 == n_2);
        swap(labels[n_1], labels[n_2]);
    }
    return;
}

//Selection operator
int Genetics::selection_operator(vector<vector<int>> &population, Random &rnd, int p) {
    double r = rnd.Rannyu();
    int j = int(population.size() * pow(r, p));
    return j;
}

//Shift operator
void Genetics::shift_operator(double prob, vector<int> &labels, int N_elem, int shift, Random &rnd) {
    if (rnd.Rannyu() < prob) {
        vector<int> labels_elem(N_elem, 0);
        vector<int> last_labels(shift, 0);
        int n_1 = int(rnd.Rannyu(2., 70.));
        for (int i = 0; i < N_elem; i++) {
            labels_elem[i] = labels[n_1 + i];
        }
        for (int i = 0; i < shift; i++) {
            last_labels[i] = labels[n_1 + N_elem + i];
        }
        for (int i = 0; i < shift; i++) {
            labels[n_1 + i] = last_labels[i];
        }
        for (int i = 0; i < N_elem; i++) {
            labels[n_1 + shift + i] = labels_elem[i];
        }
    }
    return;
}

//M-permutation
void Genetics::m_permutation(double prob, vector<int> &labels, Random &rnd) {
    if (rnd.Rannyu() < prob) {
        int m = rnd.Rannyu(2, n_cities / 2);
        int start1 = rnd.Rannyu(1, n_cities - 2 * m);
        int start2 = rnd.Rannyu(1, n_cities - m);
        while (abs(start1 - start2) < m) {
            start2 = rnd.Rannyu(1, n_cities - m);
        }
        for (int i = 0; i < m; ++i) {
            swap(labels[start1 + i], labels[start2 + i]);
        }
    }
    return;
}

//Inverse operator
void Genetics::inverse_operator(double prob, vector<int> &labels, Random &rnd) {
    if (rnd.Rannyu() < prob) {
        int m = int(rnd.Rannyu(2, n_cities));
        int start = int(rnd.Rannyu(1, n_cities - m));
        reverse(labels.begin() + start, labels.begin() + start + m);
    }
    return;
}

pair<vector<int>, vector<int>>
Genetics::cross_over_operator(vector<int> &parent_1, vector<int> &parent_2, Random &rnd) {
    int len = parent_1.size();

    // Scegli un punto di crossover casuale, escludendo l'ultimo elemento
    int crossover_point = int(rnd.Rannyu(1., len - 1));

    // Array per i figli, escluso l'ultimo elemento
    std::vector<int> offspring1(parent_1.begin(), parent_1.begin() + crossover_point);
    std::vector<int> offspring2(parent_2.begin(), parent_2.begin() + crossover_point);
    // Crea un insieme per tenere traccia degli elementi già presenti
    std::unordered_set<int> offspring1_set(offspring1.begin(), offspring1.end());
    std::unordered_set<int> offspring2_set(offspring2.begin(), offspring2.end());

    // Aggiungi gli elementi mancanti da parent2 a offspring1 mantenendo l'ordine, escludendo l'ultimo elemento
    for (int i = 0; i < len - 1; ++i) {
        int elem = parent_2[i];
        if (offspring1_set.find(elem) == offspring1_set.end()) {
            offspring1.push_back(elem);
            offspring1_set.insert(elem);
        }
    }

    // Aggiungi gli elementi mancanti da parent1 a offspring2 mantenendo l'ordine, escludendo l'ultimo elemento
    for (int i = 0; i < len - 1; ++i) {
        int elem = parent_1[i];
        if (offspring2_set.find(elem) == offspring2_set.end()) {
            offspring2.push_back(elem);
            offspring2_set.insert(elem);
        }
    }

    // Aggiungi l'ultimo elemento invariato ai figli
    offspring1.push_back(parent_1.back());
    offspring2.push_back(parent_2.back());

    return {offspring1, offspring2};
}

//Computing distances
double Genetics::compute_best_path(vector<int> &labels) {
    double path_length = 0.;
    for (int j = 0; j <= n_cities - 1; j++) {
        path_length += path.L1_norm(path.getCity(labels[j]), path.getCity(labels[j + 1]));
    }
    return path_length;
}

//Computing average of the best half of the population
double Genetics::compute_half_best_path(vector<vector<int>> &population) {
    double sum_l = 0.;
    for (int i = 0; i < pop_size / 2; i++) {
        sum_l += compute_best_path(population[i]);
    }
    return sum_l / (pop_size / 2);
}

//Get coordinates of the city
vector<double> Genetics::getProvinceCoordinates(int i) {
    return path.getCity(i).getProvinceCoordinates();
}


