/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include "functions.h"
#include <vector>
#include "Path.h"
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    field<Path> population;
    const int pop_size = 100;
    const int n_cities = 34;
    vector<double> coordinates(n_cities, 0.);
    vector<int> labels(n_cities + 1, 0);
    double r = 1.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    /*
     CREATION OF THE POPULATION
    */
    population.set_size(pop_size);

    //GENERATION OF THE LABELS AND COORDINATES
    labels[0] = 1;
    labels[n_cities] = 1;
    for (int k = 1; k < n_cities; k++) {
        labels[k] = k + 1;
    }
    for (int j = 0; j < n_cities; j++) {
        coordinates[j] = rnd.Rannyu(0., 2 * M_PI);
    }
    //GENERATION OF THE POPULATION
    for (int i = 0; i < pop_size; i++) {
        population(i).initialize(n_cities + 1);
        population(i).setCity(labels[0], coordinates[0], 0);
        population(i).setCity(labels[0], coordinates[0], n_cities);
        for (int j = 1; j < n_cities; j++) {
            population(i).setCity(labels[j], coordinates[j], j);
        }
    }

    /*
     * PAIR MUTATION OPERATION: FIRST OPERATION
     */
    int count = 1;
    for (int i = 1; i <= (n_cities - 1) / 2; i++) {
        int label_i = population(i).getCity(count).getLabel();
        int label_j = population(i).getCity(count + 1).getLabel();
        double coordinate_i = population(i).getCity(count).getCoordinate();
        double coordinate_j = population(i).getCity(count + 1).getCoordinate();
        swap(label_i, label_j);
        swap(coordinate_i, coordinate_j);
        population(i).setCity(label_i, coordinate_i, count);
        population(i).setCity(label_j, coordinate_j, count + 1);
        count += 2;
    }

    /*
     * PAIR MUTATION OPERATION: SECOND OPERATION ----> RANDOM PAIR MUTATION
     */

    for (int i = 0; i < population.size(); i++) {
        for (int n = 1; n <= pop_size; n++) {
            int n_1 = (rand() % (33 - 1 + 1)) + 1;
            int n_2 = (rand() % (33 - 1 + 1)) + 1;
            int label_i = population(i).getCity(n_1).getLabel();
            int label_j = population(i).getCity(n_2).getLabel();
            double coordinate_i = population(i).getCity(n_1).getCoordinate();
            double coordinate_j = population(i).getCity(n_2).getCoordinate();
            swap(label_i, label_j);
            swap(coordinate_i, coordinate_j);
            population(i).setCity(label_i, coordinate_i, n_1);
            population(i).setCity(label_j, coordinate_j, n_2);
        }
    }

    //THE STARTING POPULATION IS READY

    for (int i = 0; i < population.size(); i++) {
        for (int k = 0; k <= n_cities; k++) {
            cout << population(i).getCity(k).getLabel() << " ";
        }
        cout << endl;
    }

    //CHECK THE STARTING POPULATION FULFILS THE BONDS
    for (int i = 0; i < pop_size; i++) {
        vector<int> u = population(i).getLabels();
        //cout << check_function(u) << endl;
    }

    //SETTING PATH LENGTH FOR EACH PATH OF THE POPULATION
    for (int i = 0; i < pop_size; i++) {
        population(i).setPathLength(0.);
        double path = population(i).getPathLength();
        for (int k = 0; k <= n_cities - 1; k++) {
            path += L1_norm(population(i).getCity(k), population(i).getCity(k + 1), r);
        }
        population(i).setPathLength(path);
    }

    //ORDER THE POPULATION ACCORDING TO A FITNESS BASIS (FROM THE LONGEST PATH TO THE SHORTEST PATH)

    // Convertire il field in un vector per utilizzare std::stable_sort
    vector<Path> population_vector;
    for (size_t i = 0; i < population.n_elem; ++i) {
        population_vector.push_back(population(i, 0));
    }

    // Utilizziamo stable_sort con una lambda function
    //SISTEMARE L'ALGORITMO DI SORTING
    stable_sort(population_vector.begin(), population_vector.end(), []( Path &a,  Path &b) {
        return a.getPathLength() < b.getPathLength();
    });

    // Copia degli elementi ordinati di nuovo nel field
    for (size_t i = 0; i < population_vector.size(); ++i) {
        population(i, 0) = population_vector[i];
    }

    // Stampa delle persone ordinate
    for (size_t i = 0; i < population.n_elem; ++i) {
        cout << population(i).getPathLength() << endl;
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
