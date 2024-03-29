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
    vector<double> average(N, 0.);
    vector<double> error(N, 0.);
    tuple<vector<vector<double>>, vector<int>> discrete_case;
    tuple<vector<vector<double>>, vector<int>> continuum_case;
    //vector<vector<double>> sum_c(N, vector<double>(L, 0.));
    vector<vector<double>> ave(N, vector<double>(L, 0.));
    vector<vector<double>> ave2(N, vector<double>(L, 0.));
    //vector<double> count_c(n, 0.);
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

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < L; i++) {
            ave[k][i] = (get<0>(discrete_case)[k][i] / L);
            ave2[k][i] = pow(ave[k][i], 2);
        }
    }
    for (int i = 0; i < L; i++) {
        double sum1 = 0.;
        double sum2 = 0.;
        for (int k = 0; k < N; k++) {
            sum1 += ave2[k][i];
            sum2 += ave[k][i];
        }
        average[i] = sum2 / N;
        error[i] = sqrt(((sum1 / N) - pow(sum2 / N, 2)) / N);
    }

    writeOnFile(average, error, "results_1.dat");

    //Continuum case
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            count_c[0] = 0.;
            count_c[1] = 0.;
            count_c[2] = 0.;
            for (int i = 0; i < L; i++) {
                double theta = rnd.Rannyu(0, M_PI);
                double phi = rnd.Rannyu(0, 2 * M_PI);
                count_c[0] += sin(theta) * sin(phi);
                count_c[1] += sin(theta) * cos(phi);
                count_c[2] += cos(theta);
                sum_c[k][i] += sqrt((pow(count_c[0], 2) + pow(count_c[1], 2) + pow(count_c[2], 2)));
            }
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < L; i++) {
            ave[k][i] = (sum_c[k][i] / L);
            ave2[k][i] = pow(ave[k][i], 2);
        }
    }
    for (int i = 0; i < L; i++) {
        double sum1 = 0.;
        double sum2 = 0.;
        for (int k = 0; k < N; k++) {
            sum1 += ave2[k][i];
            sum2 += ave[k][i];
        }
        average[i] = sum2 / N;
        error[i] = sqrt(((sum1 / N) - pow(sum2 / N, 2)) / N);
    }

    writeOnFile(average, error, "results_2.dat");

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
