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
#include "random.h"
#include <vector>
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
    Genetics genetics_circle;
    Genetics genetics_square;
    int pop_size = 1000;
    int n_cities = 34;
    int n_generations = 200;

    vector<int> best_path(n_cities + 1, 0);
    vector<int> father(n_cities + 1, 0);
    vector<int> mother(n_cities + 1, 0);
    pair<vector<int>, vector<int>> sons;
    vector<double> probabilities = {0.025, 0.025, 0.025, 0.025, 0.9};
    vector<vector<int>> first_pop_circle(pop_size, vector<int>(n_cities + 1, 0));
    vector<vector<int>> first_pop_square(pop_size, vector<int>(n_cities + 1, 0));
    double r = 1.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");


    //STEADY STATE GENETIC ALGORITHM
    //SOSTITUIRE LA POPOLAZIONE SANTO IDDDIOOOO!!!!
    //CIRCLE
    genetics_circle.setPopSize(pop_size);
    genetics_circle.setCitiesPath(n_cities);
    genetics_circle.setProbabilities(probabilities);
    genetics_circle.initialize_path_circle(rnd);
    first_pop_circle = genetics_circle.first_pop(rnd);

    ofstream WriteResults;
    WriteResults.open("average_circle.dat");
    ofstream WriteResults1;
    WriteResults1.open("best_path_circle_coordinates.dat");

    for (int i = 0; i < n_generations; i++) {
        sons.first = father;
        sons.second = mother;
        genetics_circle.sort_paths(first_pop_circle);
        if (rnd.Rannyu() < genetics_circle.getProbabilities()[4]) {
            father = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
            mother = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
            sons = genetics_circle.cross_over_operator(father, mother, rnd);
            genetics_circle.check_function(sons.first);
            genetics_circle.check_function(sons.second);
            first_pop_circle[pop_size - 1] = sons.first;
            first_pop_circle[pop_size - 2] = sons.second;
        }
        genetics_circle.sort_paths(first_pop_circle);
        best_path = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
        genetics_circle.pair_permutation(rnd.Rannyu(), best_path, rnd);
        genetics_circle.check_function(best_path);
        first_pop_circle[pop_size - 1] = best_path;
        genetics_circle.sort_paths(first_pop_circle);
        best_path = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
        genetics_circle.shift_operator(rnd.Rannyu(), best_path, int(rnd.Rannyu(1, 4)), int(rnd.Rannyu(1, 3)), rnd);
        genetics_circle.check_function(best_path);
        first_pop_circle[pop_size - 1] = best_path;
        genetics_circle.sort_paths(first_pop_circle);
        best_path = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
        genetics_circle.m_permutation(rnd.Rannyu(), best_path, int(rnd.Rannyu(1, 6)), rnd);
        genetics_circle.check_function(best_path);
        first_pop_circle[pop_size - 1] = best_path;
        genetics_circle.sort_paths(first_pop_circle);
        best_path = genetics_circle.selection_operator(first_pop_circle, rnd, 3);
        genetics_circle.inverse_operator(rnd.Rannyu(), best_path, int(rnd.Rannyu(1, 6)), rnd);
        genetics_circle.check_function(best_path);
        first_pop_circle[pop_size - 1] = best_path;
        genetics_circle.sort_paths(first_pop_circle);


        if (WriteResults.is_open()) {
            WriteResults << i << " " << genetics_circle.compute_best_path(first_pop_circle[0], r) << " "
                         << genetics_circle.compute_half_best_path(first_pop_circle, r) << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;

    }


    for (int i = 0; i < pop_size; i++) {
        for (int k = 0; k <= n_cities; k++) {
            cout << first_pop_circle[i][k] << " ";
        }
        cout << endl;
    }

    //Printing the coordinates of the best path in (x,y) cartesian coordinates
    for (int k = 0; k <= n_cities - 1; k++) {
        if (WriteResults1.is_open()) {
            WriteResults1 << r * cos(genetics_circle.getCityCoordinate(first_pop_circle[0][k])) << " "
                          << r * sin(genetics_circle.getCityCoordinate(first_pop_circle[0][k])) << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
    }

    /*
    //SQUARE
    genetics_square.setPopSize(pop_size);
    genetics_square.setCitiesPath(n_cities);
    genetics_square.setProbabilities(probabilities);
    genetics_square.initialize_path_square(rnd);
    first_pop_square = genetics_square.first_pop(rnd);

    ofstream WriteResults2;
    WriteResults2.open("average_square.dat");
    ofstream WriteResults3;
    WriteResults3.open("best_path_square_coordinates.dat");


    for (int i = 0; i < n_generations; i++) {
        best.first = best_path;
        best.second = j;
        father_pair.first = father;
        father_pair.second = j_father;
        mother_pair.first = mother;
        mother_pair.second = j_mother;
        sons.first = father;
        sons.second = mother;
        genetics_square.sort_paths(first_pop_square);
         if (rnd.Rannyu() < genetics_square.getProbabilities()[4]) {
            father_pair = genetics_square.selection_operator(first_pop_square, rnd, 3);
            mother_pair = genetics_square.selection_operator(first_pop_square, rnd, 3);
            sons = genetics_square.cross_over_operator(father_pair.first, mother_pair.first);
            genetics_square.check_function(sons.first);
            genetics_square.check_function(sons.second);
            first_pop_square[pop_size - 1] = sons.first;
            first_pop_square[pop_size - 2] = sons.second;
        }
        genetics_square.sort_paths(first_pop_square);
        best = genetics_square.selection_operator(first_pop_square, rnd, 3);
        genetics_square.pair_permutation(rnd.Rannyu(), best.first);
        genetics_square.shift_operator(rnd.Rannyu(), best.first, (rand() % (3 - 1 + 1)) + 1,
                                       (rand() % (2 - 1 + 1)) + 1);
        genetics_square.m_permutation(rnd.Rannyu(), best.first, (rand() % (5 - 1 + 1)) + 1);
        genetics_square.inverse_operator(rnd.Rannyu(), best.first, (rand() % (5 - 1 + 1)) + 1);
        genetics_square.check_function(best.first);
        first_pop_square[pop_size - 1] = best.first;
        genetics_square.sort_paths(first_pop_square);




        if (WriteResults2.is_open()) {
            WriteResults2 << i << " " << genetics_square.compute_best_path(first_pop_circle[0], r) << " "
                         << genetics_circle.compute_half_best_path(first_pop_circle, r) << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;

    }

    //Printing the coordinates of the best path in (x,y) cartesian coordinates
    for (int k = 0; k <= n_cities - 1; k++) {
        if (WriteResults3.is_open()) {
            WriteResults3 << << " "
                          << << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
    }

    */
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
