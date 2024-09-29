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
    int M = 10000; //Throws
    int SA_steps = 1000;
    int N = 100; //Blocks
    double x = 1.; //Starting point
    double metropolis_step = 2.8;
    double mu = 4.;
    double sigma = 4.;
    double mu_best = 0.;
    double sigma_best = 0.;
    double H_min = 0.;
    double beta = 1.;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    tuple<vector<double>, vector<double>> metropolis_old;
    tuple<vector<double>, vector<double>> metropolis_new;
    tuple<vector<double>, vector<double>> H_old;
    tuple<vector<double>, vector<double>> H_new;
    tuple<vector<double>, vector<double>> hamiltonian_GS;

    ofstream WriteResults;
    WriteResults.open("results_H_vs_SA_steps.dat");
    for (int i = 0; i < SA_steps; i++) {
        double deltaSA = pow(double(beta), -0.5);
        metropolis_old = Metropolis_Uniform(x, rnd, metropolis_step, psi_trial, potential, kinetic_energy,
                                            mu, sigma,
                                            N,
                                            M);
        H_old = cumulativeAverage(get<0>(metropolis_old), get<1>(metropolis_old));
        double mu_new = rnd.Rannyu(mu - deltaSA, mu + deltaSA);
        double sigma_new = rnd.Rannyu(sigma - deltaSA, sigma + deltaSA);
        metropolis_new = Metropolis_Uniform(x, rnd, metropolis_step, psi_trial, potential, kinetic_energy,
                                            mu_new, sigma_new,
                                            N,
                                            M);
        H_new = cumulativeAverage(get<0>(metropolis_new), get<1>(metropolis_new));
        double p = 1.0;
        if (get<0>(H_new)[N - 1] > get<0>(H_old)[N - 1])
            p = exp((-1.) * beta * (get<0>(H_new)[N - 1] - get<0>(H_old)[N - 1]));
        double rand = rnd.Rannyu();
        if (rand <= p) {
            mu = mu_new;
            sigma = sigma_new;
            H_old = H_new;
        }
        if (WriteResults.is_open()) {
            WriteResults << get<0>(H_old)[N - 1] << " " << get<1>(H_old)[N - 1] << " " << mu << " " << sigma << " "
                         << beta << " " << "\t" << endl;
        } else cerr << "PROBLEM: Unable to open random.out" << endl;
        beta += 0.5;
        if (get<0>(H_old)[N - 1] < H_min) {
            mu_best = mu;
            sigma_best = sigma;
            H_min = get<0>(H_old)[N - 1];
        }
    }

    cout << mu_best << " " << sigma_best << " " << H_min << endl;

    //Evaluation of the ground state energy with the optimal parameters
    hamiltonian_GS = Metropolis_Uniform(x, rnd, metropolis_step, psi_trial, potential, kinetic_energy, mu_best,
                                        sigma_best, N,
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
