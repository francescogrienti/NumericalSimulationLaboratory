/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Genetics.h"
#include <fstream>
#include <iostream>
#include "Path.h"
#include "random.h"
#include <vector>
#include <algorithm>
#include <cmath>

//RANDOM NUMBER GENERATOR INITIALIZATION

Random initialize(Random rnd, vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file) {
    ifstream Primes(prime_file);
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input(input_file);
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    return rnd;
}

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    Genetics genetics;
    int pop_size = 100;
    int n_cities = 34;
    int n_generations = 500000;

    vector<int> best_path(n_cities + 1, 0);
    int j = 0;
    pair<vector<int>, int> best;
    vector<int> father(n_cities + 1, 0);
    int j_father = 0;
    pair<vector<int>, int> father_pair;
    vector<int> mother(n_cities + 1, 0);
    int j_mother = 0;
    pair<vector<int>, int> mother_pair;
    pair<vector<int>, vector<int>> sons;
    vector<double> probabilities = {0.1, 0.1, 0.1, 0.1, 0.6};
    vector<vector<int>> first_pop(pop_size, vector<int>(n_cities + 1, 0));
    double r = 1.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");
    genetics.setPopSize(pop_size);
    genetics.setCitiesPath(n_cities);
    genetics.setProbabilities(probabilities);
    genetics.initialize_path(rnd);
    first_pop = genetics.first_pop(rnd);

    ofstream WriteResults;
    WriteResults.open("average_circle.dat");
    ofstream WriteResults1;
    WriteResults1.open("best_path_circle_coordinates.dat");

    for (int i = 0; i < n_generations; i++) {
        best.first = best_path;
        best.second = j;
        father_pair.first = father;
        father_pair.second = j_father;
        mother_pair.first = mother;
        mother_pair.second = j_mother;
        sons.first = father;
        sons.second = mother;
        genetics.sort_paths(first_pop);
        best = genetics.selection_operator(first_pop, rnd, 3);
        genetics.pair_permutation(rnd.Rannyu(), best.first);
        genetics.shift_operator(rnd.Rannyu(), best.first, (rand() % (3 - 1 + 1)) + 1, (rand() % (2 - 1 + 1)) + 1);
        genetics.m_permutation(rnd.Rannyu(), best.first, (rand() % (5 - 1 + 1)) + 1);
        genetics.inverse_operator(rnd.Rannyu(), best.first, (rand() % (5 - 1 + 1)) + 1);
        genetics.check_function(best.first);
        first_pop[pop_size - 1] = best.first;
        genetics.sort_paths(first_pop);
        if (rnd.Rannyu() < genetics.getProbabilities()[4]) {
            father_pair = genetics.selection_operator(first_pop, rnd, 3);
            mother_pair = genetics.selection_operator(first_pop, rnd, 3);
            sons = genetics.cross_over_operator(father_pair.first, mother_pair.first);
            genetics.check_function(sons.first);
            genetics.check_function(sons.second);
            first_pop[pop_size - 1] = sons.first;
            first_pop[pop_size - 2] = sons.second;
        }
        genetics.sort_paths(first_pop);

        /*
        if (WriteResults.is_open()) {
            WriteResults << i << " " << genetics.compute_best_path(first_pop[0], r) << " "
                         << genetics.compute_half_best_path(first_pop, r) << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
        */
    }

    //Printing the coordinates of the best path in (x,y) cartesian coordinates
    for (int k = 0; k <= n_cities - 1; k++) {
        if (WriteResults1.is_open()) {
            WriteResults1 << r * cos(genetics.getCityCoordinate(first_pop[0][k])) << " "
                          << r * sin(genetics.getCityCoordinate(first_pop[0][k])) << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
    }


    rnd.SaveSeed();


    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
