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
#include "City.h"
#include "Population.h"

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    Population population;
    Path path;
    const int pop_size = 100;
    const int n_cities = 34;
    vector<double> coordinates(n_cities, 0.);
    vector<int> labels(n_cities, 0);
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //CREATION OF THE POPULATION

    //CREATION OF THE FIRST OBJECT OF TYPE PATH
    labels[0] = 1;
    labels[n_cities-1] = 1;
    for (int k = 2; k < n_cities-1; k++) {
        labels[k] = k;
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
