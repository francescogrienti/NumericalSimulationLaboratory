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
#include <algorithm>


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
            //cout << p1 << " " << p2 << endl;
        } else if (rank == 1) {
            int linesToSkip = 1;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
            //cout << p1 << " " << p2 << endl;
        } else if (rank == 2) {
            int linesToSkip = 2;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
            //cout << p1 << " " << p2 << endl;
        } else if (rank == 3) {
            int linesToSkip = 3;  // Numero di righe da saltare
            skipLines(Primes, linesToSkip);
            Primes >> p1 >> p2;
            //cout << p1 << " " << p2 << endl;
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

    //MPI Initialization
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;

    //Initializing variables
    Random rnd;
    Genetics evolution;
    int pop_size = 1200;
    int n_cities = 110;
    int n_generations = 1000;
    vector<int> father(n_cities + 1, 0);
    vector<int> mother(n_cities + 1, 0);
    pair<vector<int>, vector<int>> sons;
    vector<vector<int>> first_population(pop_size, vector<int>(n_cities + 1, 0));
    vector<vector<int>> evo_population(pop_size, vector<int>(n_cities + 1, 0));
    vector<double> probab = {0.07, 0.07, 0.07, 0.07, 0.9};
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in", rank);

    //Variables for migration
    int N_migr = 3;
    vector<int> message(n_cities + 1, 0);
    vector<int> message2(n_cities + 1);  // vettore d'appoggio per la migrazione
    int itag1 = 1;
    int itag2 = 2;
    vector<int> ranks = {0, 1, 2, 3};


    //Setting
    evolution.setPopSize(pop_size);
    evolution.setCitiesPath(n_cities);
    evolution.initialize_path("cap_prov_ita.dat");
    first_population = evolution.first_pop(rnd);

    ofstream WriteResults;
    WriteResults.open("average_" + to_string(rank) + ".dat");
    ofstream WriteResults1;
    WriteResults1.open("best_path_" + to_string(rank) + "_coordinates.dat");

    //Running of the genetic algorithm
    for (int i = 0; i < n_generations; i++) {
        evolution.sort_paths(first_population);
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
                evolution.shift_operator(probab[3], sons.first, int(rnd.Rannyu(1, 10)), int(rnd.Rannyu(1, 10)),
                                         rnd);
                evolution.shift_operator(probab[3], sons.second, int(rnd.Rannyu(1, 10)), int(rnd.Rannyu(1, 10)),
                                         rnd);
                evo_population[2 * k] = sons.first;
                evo_population[2 * k + 1] = sons.second;
            } else {
                evo_population[2 * k] = father;
                evo_population[2 * k + 1] = mother;
            }
        }
        evolution.sort_paths(evo_population);
        //MIGRATIONS
        if (i % N_migr == 0) {
            random_shuffle(ranks.begin(), ranks.end());
            for (int j = 0; j <= n_cities; j++) {
                message[j] = evo_population[0][j]; //Best individual
                message2[j] = evo_population[0][j]; // Prendo il miglior individuo della popolazione locale
            }
            if (rank == ranks[1]) {
                MPI_Send(&message[0], n_cities + 1, MPI_INTEGER, ranks[0], itag1,
                         MPI_COMM_WORLD);
                MPI_Recv(&message2[0], n_cities + 1, MPI_INTEGER, ranks[0], itag2, MPI_COMM_WORLD,
                         &stat);
            } else if (rank == ranks[0]) {
                MPI_Send(&message2[0], n_cities + 1, MPI_INTEGER, ranks[1], itag2,
                         MPI_COMM_WORLD);
                MPI_Recv(&message[0], n_cities + 1, MPI_INTEGER, ranks[1], itag1, MPI_COMM_WORLD,
                         &stat);
            }
            // Sostituisco gli individui locali con quelli ricevuti
            if (rank == ranks[1]) {
                for (int j = 0; j <= n_cities; j++) evo_population[0][j] = message2[j];
            } else if (rank == ranks[0]) {
                for (int j = 0; j <= n_cities; j++) evo_population[0][j] = message[j];
            }
            evolution.sort_paths(evo_population);
        }

        first_population = evo_population;
        if (WriteResults.is_open()) {
            WriteResults << i << " " << evolution.compute_best_path(first_population[0]) << " "
                         << evolution.compute_half_best_path(first_population) << " " << "\t"
                         << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;

        //Printing the coordinates of the best path in (x,y) cartesian coordinates
        for (int k = 0; k <= n_cities; k++) {
            if (WriteResults1.is_open()) {
                WriteResults1 << evolution.getProvinceCoordinates(first_population[0][k])[0] << " "
                              << evolution.getProvinceCoordinates(first_population[0][k])[1] << " " << "\t"
                              << endl;
            } else cerr << "PROBLEM: Unable to open random.out" << endl;
        }

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
