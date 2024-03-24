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
#include "vector"

using namespace std;

int main(int argc, char *argv[]) {

    //Variables declaration
    Random rnd;
    const int M = 1000000;
    const vector<int> N = {1, 2, 10, 100};
    const double lambda = 1.;
    const double mu = 0.;
    const double gamma = 1.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");
    writeOnFile(M, N, rnd, lambda, expon_prob_sum, "results_1.dat");
    writeOnFile(M, N, rnd, mu, gamma, cauchy_prob_sum, "results_2.dat");
    writeOnFile(M, N, rnd, uniform_prob_sum, "results_3.dat");

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
