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
#include <vector>
#include <tuple>
#include "random.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]) {

    //Variables declaration
    Random rnd;
    const int M = 1000000; //Number of total throws
    const int N = 100; //Number of blocks
    const double mu = 0.5; //Expected mean
    const int expecValue = M / N;
    tuple<vector<double>, vector<double>> averages;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives;
    vector<double> chi2(N, 0.);
    vector<int> count(N, 0);
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //Evaluation of the mean
    averages = mean(M, N, rnd);
    cumulatives = cumulativeAverage(get<0>(averages), get<1>(averages));
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_1.dat");

    //Evaluation of the standard deviation
    averages = mean(M, N, rnd, mu);
    cumulatives = cumulativeAverage(get<0>(averages), get<1>(averages));
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_2.dat");

    //Chi-2
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < M; k++) {
            double r = rnd.Rannyu();
            count[(int) (r * N)] += 1;
        }
        for (int j = 0; j < N; j++) {
            chi2[j] += double(pow(count[j] - expecValue, 2)) / (double) (expecValue);
            count[j] = 0;
        }
    }

    writeOnFile(chi2, "results_3.dat");

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
