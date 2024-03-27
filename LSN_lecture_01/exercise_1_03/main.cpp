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
    const int M = 1000000; //Throws
    const int N = 100; //Blocks
    const int L = M / N; //Throws per block
    const double l = 1.; //Length of the needle
    const double d = 2.; //Spacing among lines
    int N_hit = 0;
    vector<double> ave(N, 0.);
    vector<double> ave2(N, 0.);
    tuple<vector<double>, vector<double>, vector<double >> cumulatives;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    for (int i = 0; i < N; i++) {
        N_hit = 0;
        for (int j = 0; j < L; j++) {
            double x1 = rnd.Rannyu(0, d);
            double x2 = 0.;
            double angle = rnd.Theta();
            x2 = x1 + l * cos(2 * angle);
            if (x2 < 0 || x2 > d) {
                N_hit += 1;
            }
            ave[i] = (2.0 * l * L) / (N_hit * d); //Store average values for each block
            ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
        }
    }

    cumulatives = cumulativeAverage(ave, ave2);
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_1.dat");
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
