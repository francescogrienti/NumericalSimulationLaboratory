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
    int SA_steps = 10000; //Throws
    int N = 200; //Blocks
    double x = 1.1; //Starting point
    double x_new = 0.;
    const int n_iter = 5;
    double metropolis_step = 2.9;
    double mu = 1.;
    double sigma = 0.;
    double T_start = 0.1;
    double delta_mu = 0.1;
    double delta_sigma = 0.1;
    double delta_T = 0.0005;
    double x_delta = 2.5;
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

    for (double T = T_start; T >= delta_T; T -= delta_T) {
        tuple<vector<double>, vector<double>> metropolis_old;
        tuple<vector<double>, vector<double>> metropolis_new;
        tuple<vector<double>, vector<double>> H_old;
        tuple<vector<double>, vector<double>> H_new;
        double beta = 1. / T;
        metropolis_step = rnd.Rannyu(metropolis_step - 1., metropolis_step + 1.);
        for (int i = 0; i < n_iter; i++) {
            metropolis_old = Metropolis_Uniform(x, rnd, metropolis_step, pdf_function, potential, kinetic_energy,
                                                mu, sigma,
                                                N,
                                                SA_steps, "results_x_hamiltonian_GS.dat");
            H_old = cumulativeAverage(get<0>(metropolis_old), get<1>(metropolis_old), "results_energy_GS.dat");
            x_new = rnd.Rannyu(x - x_delta, x + x_delta);
            mu_new = fabs(mu + delta_mu * (rnd.Rannyu() - 0.5));
            sigma_new = fabs(sigma + delta_sigma * (rnd.Rannyu() - 0.5));
            metropolis_new = Metropolis_Uniform(x_new, rnd, metropolis_step, pdf_function, potential, kinetic_energy,
                                                mu_new, sigma_new,
                                                N,
                                                SA_steps, "results_x_hamiltonian_GS.dat");
            H_new = cumulativeAverage(get<0>(metropolis_new), get<1>(metropolis_new), "results_energy_GS.dat");
            double p = 1.0;
            if (get<0>(H_new)[N - 1] > get<0>(H_old)[N - 1])
                p = exp(-beta * (get<0>(H_new)[N - 1] - get<0>(H_old)[N - 1]));
            if (p >= rnd.Rannyu()) {
                mu = mu_new;
                sigma = sigma_new;
                get<0>(H_old) = get<0>(H_new);
            }
        }
        if (WriteResults.is_open()) {
            WriteResults << get<0>(H_old)[N - 1] << " " << get<1>(H_old)[N - 1] << " " << "\t" << endl;
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
