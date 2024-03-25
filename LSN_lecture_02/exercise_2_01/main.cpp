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
#include <tuple>

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    const int M = 10000;
    const int N = 100;
    const int L = M / N;
    vector<double> mean_integral(N, 0.);
    vector<double> mean_integral_2(N, 0.);
    vector<double> ave(N, 0.); //Arrays for storing the average and the square of the average for each block
    vector<double> ave2(N, 0.); //Arrays for storing the average and the square of the average for each block
    tuple<vector<double>, vector<double>, vector<double>> cumulatives;
    vector<int> seed(4, 0);
    int p1 = 0;
    int p2 = 0;

    rnd = initialize(rnd, seed, p1, p2, "Primes", "seed.in");

    //Method 1 - Generate 10^4 integrals
    //Sampling a uniform probability distribution in [0,1]
    for (int i = 0; i < N; i++) {
        double sum1 = 0.;
        for (int j = 0; j < L; j++) {
            double sum2 = 0.;
            for (int k = 0; k < N; k++) {
                double r = rnd.Rannyu();
                sum2 += func_f(r);
            }
            mean_integral[j] = sum2 / N;
            sum1 += mean_integral[j];
            mean_integral[j] = 0.;
        }
        ave[i] = double(sum1 / L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
    }

    cumulatives = cumulativeAverage(ave, ave2);
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_1.dat");

    // Importance sampling - sampling a non-uniform probability distribution in [0,1]
    // using the accept-reject method for sampling the distribution
    for (int i = 0; i < N; i++) {
        double sum1 = 0.;
        for (int j = 0; j < L; j++) {
            double sum2 = 0.;
            for (int k = 0; k < N; k++) {
                double r = 0.;
                double t = 0.;
                do {
                    r = rnd.Rannyu();
                    t = rnd.Rannyu(0, 1);
                } while (r >= prob_p(t) / (1.5));
                sum2 += func_g(t);
            }
            mean_integral[j] = sum2 / N;
            sum1 += mean_integral[j];
            mean_integral[j] = 0.;
        }
        ave[i] = double(sum1 / L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
    }


    cumulatives = cumulativeAverage(ave, ave2);
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_2.dat");


    //Method 2 - Generate 10^4 points for the evaluation of 10^2 integrals
    //Sampling a uniform probability distribution in [0,1]
    for (int j = 0; j < L; j++) {
        double sum = 0.;
        for (int k = 0; k < N; k++) {
            double r = rnd.Rannyu();
            sum += func_f(r);
        }
        mean_integral[j] = sum / N;
        mean_integral_2[j] = pow(mean_integral[j], 2);
    }

    cumulatives = cumulativeAverage(mean_integral, mean_integral_2);
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_3.dat");

    // Importance sampling - sampling a non-uniform probability distribution in [0,1]
    // using the accept-reject method for sampling the distribution
    for (int j = 0; j < L; j++) {
        double sum = 0.;
        for (int k = 0; k < N; k++) {
            double r = 0.;
            double t = 0.;
            do {
                r = rnd.Rannyu();
                t = rnd.Rannyu(0, 1);
            } while (r >= prob_p(t) / (1.5));
            sum += func_g(t);
        }
        mean_integral[j] = sum / N;
        mean_integral_2[j] = pow(mean_integral[j], 2);
    }

    cumulatives = cumulativeAverage(mean_integral, mean_integral_2);
    writeOnFile(get<0>(cumulatives), get<2>(cumulatives), "results_4.dat");

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
