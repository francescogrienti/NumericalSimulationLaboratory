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
#include <tuple>

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    int M = 1000000; //Throws
    int N = 100; //Blocks
    double x = 1.; //Starting point
    tuple<vector<double>, vector<double>> hamiltonian_GS;
    double metropolis_step = 2.8;
    double mu = 0.83;
    double sigma = 0.65;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");


    hamiltonian_GS = Metropolis_Uniform(x, rnd, metropolis_step, psi_trial, potential, kinetic_energy, mu, sigma, N,
                                        M, "results_x_hamiltonian_GS.dat");

    cumulativeAverage(get<0>(hamiltonian_GS), get<1>(hamiltonian_GS), "results_energy_GS.dat");
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
