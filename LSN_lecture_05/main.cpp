/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include "random.h"
#include "functions.h"
#include <vector>
#include <tuple>

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    int M = 1000000; //Throws
    int N = 100; //Blocks
    vector<double> x_0_gauss_GS(3, 1.0);
    vector<double> x_0_gauss_ES(3, 3.0);
    vector<double> x_0_uniform_GS(3, 1.0);
    vector<double> x_0_uniform_ES(3, 3.0);
    tuple<vector<double>, vector<double>> Gauss_GS;
    tuple<vector<double>, vector<double>> Gauss_ES;
    tuple<vector<double>, vector<double>> Uniform_GS;
    tuple<vector<double>, vector<double>> Uniform_ES;
    double metropolis_step_gauss_GS = 0.75;
    double metropolis_step_gauss_ES = 1.75;
    double metropolis_step_uniform_GS = 1.2;
    double metropolis_step_uniform_ES = 1.6;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;
    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");


    Gauss_GS = Metropolis_Gauss(x_0_gauss_GS, rnd, metropolis_step_gauss_GS, pdf_wave_function_GS, N, M,
                                "results_xyz_gauss_GS.dat");
    Gauss_ES = Metropolis_Gauss(x_0_gauss_ES, rnd, metropolis_step_gauss_ES, pdf_wave_function_ES, N, M,
                                "results_xyz_gauss_ES.dat");
    Uniform_GS = Metropolis_Uniform(x_0_uniform_GS, rnd, metropolis_step_uniform_GS, pdf_wave_function_GS, N, M,
                                  "results_xyz_uniform_GS.dat");
    Uniform_ES = Metropolis_Uniform(x_0_uniform_ES, rnd, metropolis_step_uniform_ES, pdf_wave_function_ES, N, M,
                                  "results_xyz_uniform_ES.dat");

    cumulativeAverage(get<0>(Gauss_GS), get<1>(Gauss_GS), "results_radius_gauss_GS.dat");
    cumulativeAverage(get<0>(Gauss_ES), get<1>(Gauss_ES), "results_radius_gauss_ES.dat");
    cumulativeAverage(get<0>(Uniform_GS), get<1>(Uniform_GS), "results_radius_uniform_GS.dat");
    cumulativeAverage(get<0>(Uniform_ES), get<1>(Uniform_ES), "results_radius_uniform_ES.dat");


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
