//
// Created by francesco on 07/06/24.
//

#include "Genetics.h"
#include <vector>
#include <tuple>
#include <cmath>
#include "random.h"
#include <algorithm>
#include <ctime>
#include <cstdlib>

using namespace std;

Genetics::Genetics() {}
// Default constructor, does not perform any action

Genetics::~Genetics() {}
// Default destructor, does not perform any action

//Set the probability for the mutation operators
void Genetics::setProbabilities(std::vector<double> prob) {
    probabilities.resize(5);
    probabilities[0] = prob[0]; //Pair probability
    probabilities[1] = prob[1]; //Shift probability
    probabilities[2] = prob[2]; //M-probability
    probabilities[3] = prob[3]; //Inverse probability
    probabilities[4] = prob[4]; //Cross-over probability
}

//Get the probability for the mutation operators
vector<double> Genetics::getProbabilities() {
    return probabilities;

}

void Genetics::setCitiesPath(int n) {
    n_cities = n;
    return;
}

void Genetics::setPopSize(int n) {
    pop_size = n;
    return;
}

void Genetics::initialize_path(Random &rnd) {
    vector<double> coordinates(n_cities, 0.);
    vector<int> labels(n_cities, 0);
    for (int k = 0; k < n_cities; k++) {
        labels[k] = k + 1;
        coordinates[k] = rnd.Rannyu(0., 2 * M_PI);
    }
    path.initialize_path(n_cities);
    for (int k = 0; k < n_cities; k++) {
        path.setCity(labels[k], coordinates[k], k);
    }

    return;
}

//Creation of the first population using the mutation operator
vector<vector<int>> Genetics::first_pop(Random &rnd) {
    const int N_mut = 100;
    vector<vector<int>> population(pop_size, vector<int>(n_cities + 1, 0));

    //GENERATION OF THE POPULATION
    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j <= n_cities; j++) {
            population[i][0] = 1;
            population[i][n_cities] = 1;
            population[i][j + 1] = j + 2;
        }
    }
    //FIRST MUTATION: PAIR MUTATION!
    for (int i = 0; i < pop_size; i++) {
        for (int j = 1; j < N_mut; j++) {
            pair_permutation(rnd.Rannyu(), population[i]);
        }
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
void Genetics::sort_paths(std::vector<vector<int>>& population) {
    const double r = 1.0;
    vector<double> path_length(population.size(), 0.0);

    // Compute the distances of the various paths
    for (int i = 0; i < population.size(); ++i) {
        for (int j = 0; j < n_cities; ++j) {
            path_length[i] += path.L1_norm(path.getCity(population[i][j]), path.getCity(population[i][j + 1]), r);
        }
    }

    // Sorting
    // Create the struct to use for sorting
    struct VecWithVal {
        std::vector<int> vec;
        double value;
    };

    std::vector<VecWithVal> vecWithVals(population.size());
    for (int i = 0; i < population.size(); ++i) {
        vecWithVals[i].vec = population[i];
        vecWithVals[i].value = path_length[i];
    }

    auto compare = [](const VecWithVal& a, const VecWithVal& b) {
        return a.value < b.value;
    };

    std::sort(vecWithVals.begin(), vecWithVals.end(), compare);

    for (int i = 0; i < population.size(); ++i) {
        population[i] = vecWithVals[i].vec;
    }
}

void Genetics::pair_permutation(double prob, vector<int> &labels) {
    if (prob < probabilities[0]) {
        //Generation of two random integer numbers to select the pair indeces to be muted
        int n_1 = (rand() % (33 - 1 + 1)) + 1;
        int n_2 = (rand() % (33 - 1 + 1)) + 1;
        swap(labels[n_1], labels[n_2]);
    }

    return;
}

//Selection operator
pair<vector<int>, int> Genetics::selection_operator(const vector<vector<int>> &population, Random &rnd, int p) {
    double r = rnd.Rannyu();
    int j = int(population.size() * pow(r, p));
    return {population[j], j};
}

//Shift operator
//CONTROLLARE E AGGIUSTARE n_1!
void Genetics::shift_operator(double prob, vector<int> &labels, int N_elem, int shift) {
    if (prob < probabilities[0]) {
        vector<int> labels_elem(N_elem, 0);
        vector<int> last_labels(shift, 0);
        int n_1 = (rand() % (20 - 5 + 1)) + 5;
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
}

//M-permutation
//CONTROLLARE E AGGIUSTARE index1 e index2!
void Genetics::m_permutation(double prob, vector<int> &labels, int n) {
    if (prob < probabilities[2]) {
        int len = labels.size();

        if (len < 2 * n) {
            std::cerr << "Array too small to contain two subarrays of length " << n << std::endl;
            return;
        }

        std::srand(std::time(0));
        int index1 = (rand() % (20 - 5 + 1)) + 5;
        int index2;
        do {
            index2 = (rand() % (20 - 5 + 1)) + 5;
        } while (index1 == index2);

        // Scambia gli elementi dei due sotto-array
        for (int i = 0; i < n; ++i) {
            std::swap(labels[index1 + i], labels[index2 + i]);
        }
    }
}

//Inverse operator
void Genetics::inverse_operator(double prob, vector<int> &labels, int n) {
    if (prob < probabilities[3]) {
        int len = labels.size();
        int start = (rand() % (20 - 5+ 1)) + 5;

        if (len < start + n) {
            std::cerr << "Array too small to contain a subarray of length " << n << " starting at index " << start
                      << std::endl;
            return;
        }

        // Identifica l'indice di fine del sotto-array da invertire
        int end = start + n - 1;

        // Inverte il sotto-array
        while (start < end) {
            std::swap(labels[start], labels[end]);
            start++;
            end--;
        }
    }
}

pair<vector<int>, vector<int>>
Genetics::cross_over_operator(vector<int> &parent_1, vector<int> &parent_2) {
    int len = parent_1.size();

    // Inizializza la generazione di numeri casuali
    std::srand(std::time(0));

    // Scegli un punto di crossover casuale, escludendo l'ultimo elemento
    int crossover_point = std::rand() % (len - 1);

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
double Genetics::compute_best_path(vector<int> &labels, double r) {
    double path_length = 0.;
    for (int j = 0; j <= n_cities - 1; j++) {
        path_length += path.L1_norm(path.getCity(labels[j]), path.getCity(labels[j + 1]), r);
    }
    return path_length;
}

//Computing average of the best half of the population
double Genetics::compute_half_best_path(vector<vector<int>> &population, double r) {
    double sum_l = 0.;
    for (int i = 0; i < pop_size / 2; i++) {
        sum_l += compute_best_path(population[i], r);
    }
    return sum_l / (pop_size / 2);
}

//Get coordinates of the city
double Genetics::getCityCoordinate(int i) {
    return path.getCity(i).getCoordinate();
}

