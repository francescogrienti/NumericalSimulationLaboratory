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
    int N = 200; //Blocks
    double x = 1.1; //Starting point
    tuple<vector<double>, vector<double>> metro;
    tuple<vector<double>, vector<double>> hamiltonian_GS;
    double metropolis_step = 2.9;
    double mu = 1.;
    double sigma = 1.;
    double T_start = 4.;
    double T = 0.;
    double SA_steps = 350;
    double SA_mu_interval = 0.09;
    double SA_sigma_interval = 0.09;
    double T_decrem = 0.01;
    double acceptance = 0.;
    double r = 0.;
    double mu_new = 0.;
    double sigma_new = 0.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");


    ofstream WriteResults;
    WriteResults.open("results_H_vs_SA_steps.dat");
    ofstream WriteResults1;
    WriteResults1.open("results_param_vs_SA_steps.dat");
    for (int k = 0; k < SA_steps; k++) {
        metro = Metropolis_Uniform(x, rnd, metropolis_step, pdf_function, potential, kinetic_energy, mu, sigma,
                                   N,
                                   M, "results_x_hamiltonian_GS.dat");
        hamiltonian_GS = cumulativeAverage(get<0>(metro), get<1>(metro), "results_energy_GS.dat");
        if (WriteResults.is_open()) {
            WriteResults << (get<0>(hamiltonian_GS)[N - 1]) << " " << (get<1>(hamiltonian_GS)[N - 1]) << " " << "\t"
                         << endl;
        }
        if (WriteResults1.is_open()) {
            WriteResults1 << mu << " " << sigma << " " << "\t" << endl;
        }
        T = T_start - 1. * T_decrem;
        mu_new = rnd.Rannyu(mu - SA_mu_interval, x + SA_mu_interval);
        sigma_new = rnd.Rannyu(sigma - SA_sigma_interval, sigma + SA_sigma_interval);
        acceptance = min(1., exp((-1. / T) * ((potential(x) + kinetic_energy(x, mu_new, sigma_new)))) /
                             exp((-1. / T_start) * (get<0>(hamiltonian_GS)[N - 1])));
        r = rnd.Rannyu();
        if (r <= acceptance) {
            mu = mu_new;
            sigma = sigma_new;
        }
        cout << T << endl;
        T_start = T;

    }
    cout << mu << " " << sigma << endl;
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
