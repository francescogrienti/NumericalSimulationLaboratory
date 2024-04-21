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
    int M = 100000;
    int N = 200;
    int T = 1;
    int strike_price = 100;
    double r = 0.1;
    double sigma = 0.25;
    int S_0 = 100;
    double t_i = (double) T / N;

    tuple<vector<double>, vector<double>, vector<double>, vector<double>> averages_direct;
    tuple<vector<double>, vector<double>, vector<double>, vector<double>> averages_discretized;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives_call_direct;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives_call_discretized;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives_put_direct;
    tuple<vector<double>, vector<double>, vector<double>> cumulatives_put_discretized;

    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //Direct
    averages_direct = sampling(rnd, r, sigma, N, M, S_0, T, strike_price, t_i, "direct");

    //Discretized
    averages_discretized = sampling(rnd, r, sigma, N, M, S_0, T, strike_price, t_i, "discretized");

    cumulatives_call_direct = cumulativeAverage(get<0>(averages_direct), get<1>(averages_direct));
    cumulatives_call_discretized = cumulativeAverage(get<0>(averages_discretized), get<1>(averages_discretized));
    cumulatives_put_direct = cumulativeAverage(get<2>(averages_direct), get<3>(averages_direct));
    cumulatives_put_discretized = cumulativeAverage(get<2>(averages_discretized), get<3>(averages_discretized));

    writeOnFile(get<0>(cumulatives_call_direct), get<2>(cumulatives_call_direct), "results_call_direct.dat");
    writeOnFile(get<0>(cumulatives_call_discretized), get<2>(cumulatives_call_discretized),
                "results_call_discretized.dat");
    writeOnFile(get<0>(cumulatives_put_direct), get<2>(cumulatives_put_direct), "results_put_direct.dat");
    writeOnFile(get<0>(cumulatives_put_discretized), get<2>(cumulatives_put_discretized),
                "results_put_discretized.dat");

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
