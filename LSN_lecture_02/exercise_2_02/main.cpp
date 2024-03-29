/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include "random.h"
#include "functions.h"
#include <vector>
#include <tuple>

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    const int M = 10000; //Number of throws
    const int n = 3; //Possible directions (both positive and negative)
    const int N = 100; //Number of blocks
    const int L = M / N; //Number of throws per block;
    const int a = 1; //Lattice step
    tuple<vector<double>, vector<double>> results;
    tuple<vector<vector<double>>, vector<int>> discrete_case;
    tuple<vector<vector<double>>, vector<int>> continuum_case;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //Discrete case
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            get<1>(discrete_case)[0] = 0;
            get<1>(discrete_case)[1] = 0;
            get<1>(discrete_case)[2] = 0;
            for (int i = 0; i < L; i++) {
                double r = rnd.Rannyu(-1, 1);
                if (r < 0) {
                    get<1>(discrete_case)[(int) ((-1) * r * n)] -= a;
                } else if (r > 0) {
                    get<1>(discrete_case)[(int) (r * n)] += a;
                }
                get<0>(discrete_case)[k][i] += sqrt(
                        (pow(get<1>(discrete_case)[0], 2) + pow(get<1>(discrete_case)[1], 2) +
                         pow(get<1>(discrete_case)[2], 2)));
            }
        }
    }
    results = mean_and_error(get<0>(discrete_case), N, L);
    writeOnFile(get<0>(results), get<1>(results), "results_1.dat");

    //Continuum case
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            get<1>(continuum_case)[0] = 0;
            get<1>(continuum_case)[1] = 0;
            get<1>(continuum_case)[2] = 0;
            for (int i = 0; i < L; i++) {
                double theta = rnd.Rannyu(0, M_PI);
                double phi = rnd.Rannyu(0, 2 * M_PI);
                get<1>(continuum_case)[0] += sin(theta) * sin(phi);
                get<1>(continuum_case)[1] += sin(theta) * cos(phi);
                get<1>(continuum_case)[2] += cos(theta);
                get<0>(continuum_case)[k][i] += sqrt(
                        (pow(get<1>(continuum_case)[0], 2) + pow(get<1>(continuum_case)[1], 2) +
                         pow(get<1>(continuum_case)[2], 2)));
            }
        }
    }
    results = mean_and_error(get<0>(continuum_case), N, L);
    writeOnFile(get<0>(results), get<1>(results), "results_2.dat");

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
