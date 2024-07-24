/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"
#include "Genetics.h"
#include <fstream>
#include <iostream>
#include "random.h"
#include <vector>
#include <cmath>


// Funzione helper per saltare un numero specifico di righe
void skipLines(ifstream &file, int numLines) {
    string dummy;
    for (int i = 0; i < numLines; ++i) {
        if (!getline(file, dummy)) {
            cerr << "PROBLEM: Unable to skip line " << i + 1 << endl;
            break;
        }
    }
}

//RANDOM NUMBER GENERATOR INITIALIZATION
Random initialize(Random rnd, vector<int> seed, int p1, int p2, const std::string &prime_file, std::string input_file,
                  int rank) {
    ifstream Primes(prime_file);
    if (Primes.is_open()) {
        if (rank == 0) {
            Primes >> p1 >> p2;
        } else if (rank == 1) {
            int linesToSkip = 1;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
        } else if (rank == 2) {
            int linesToSkip = 2;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
        } else if (rank == 3) {
            int linesToSkip = 3;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
        }
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

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <value>" << std::endl;
        return 1; // indicate error
    }

    // argv[0] is the name of the program itself
    // argv[1] is the first argument passed by the user

    // Convert the argument to an integer
    string type = argv[1]; //Square or Circle minimization

    //MPI Initialization
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Initializing variables
    Random rnd;
    Genetics evolution;
    int pop_size = 500;
    int n_cities = 34;
    int n_generations = 100;
    vector<int> father(n_cities + 1, 0);
    vector<int> mother(n_cities + 1, 0);
    pair<vector<int>, vector<int>> sons;
    vector<vector<int>> first_population(pop_size, vector<int>(n_cities + 1, 0));
    vector<vector<int>> evo_population(pop_size, vector<int>(n_cities + 1, 0));
    double r = 1.;
    vector<double> probab = {0.08, 0.08, 0.08, 0.08, 0.9};
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in", rank);

    //Setting
    evolution.setPopSize(pop_size);
    evolution.setCitiesPath(n_cities);
    evolution.initialize_path(rnd, type);
    first_population = evolution.first_pop(rnd);


    ofstream WriteResults;
    WriteResults.open("average_" + type + ".dat");
    ofstream WriteResults1;
    WriteResults1.open("best_path_" + type + "_coordinates.dat");

    //Running of the genetic algorithm
    for (int i = 0; i < n_generations; i++) {
        evolution.sort_paths(first_population, type);
        for (int k = 0; k < pop_size / 2; k++) {
            int f = 0;
            int m = 0;
            do {
                f = evolution.selection_operator(first_population, rnd, 3);
                m = evolution.selection_operator(first_population, rnd, 3);
            } while (f == m);

            father = first_population[f];
            mother = first_population[m];
            if (rnd.Rannyu() < probab[4]) {
                sons = evolution.cross_over_operator(father, mother, rnd);
                evolution.check_function(sons.first);
                evolution.check_function(sons.second);
                evolution.pair_permutation(probab[0], sons.first, rnd);
                evolution.pair_permutation(probab[0], sons.second, rnd);
                evolution.inverse_operator(probab[1], sons.first, rnd);
                evolution.inverse_operator(probab[1], sons.second, rnd);
                evolution.m_permutation(probab[2], sons.first, rnd);
                evolution.m_permutation(probab[2], sons.second, rnd);
                evolution.shift_operator(probab[3], sons.first, int(rnd.Rannyu(1, 4)), int(rnd.Rannyu(1, 3)),
                                         rnd);
                evolution.shift_operator(probab[3], sons.second, int(rnd.Rannyu(1, 4)), int(rnd.Rannyu(1, 3)),
                                         rnd);
                evo_population[2 * k] = sons.first;
                evo_population[2 * k + 1] = sons.second;
            } else {
                evo_population[2 * k] = father;
                evo_population[2 * k + 1] = mother;
            }
        }

        evolution.sort_paths(evo_population, type);
        first_population = evo_population;
        if (WriteResults.is_open()) {
            WriteResults << i << " " << evolution.compute_best_path(first_population[0], r, type) << " "
                         << evolution.compute_half_best_path(first_population, r, type) << " " << "\t"
                         << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
    }

    //Printing the coordinates of the best path in (x,y) cartesian coordinates
    for (int k = 0; k <= n_cities; k++) {
        if (WriteResults1.is_open()) {
            if (type == "circle") {
                WriteResults1 << r * cos(evolution.getCityCoordinate(first_population[0][k])) << " "
                              << r * sin(evolution.getCityCoordinate(first_population[0][k])) << " " << "\t" << endl;
            } else if (type == "square") {
                WriteResults1 << evolution.getCitySquareCoordinates(first_population[0][k])[0] << " "
                              << evolution.getCitySquareCoordinates(first_population[0][k])[1] << " " << "\t" << endl;
            }
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
    }

    rnd.SaveSeed();
    MPI_Finalize();
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
