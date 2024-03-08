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
#include <string>
#include "random.h"

double* expon_prob_sum(Random rnd, int n, int throws, double lambda){
    double* sum = new double[throws]{0.};
    for(int i=0; i<throws; i++){
        for(int j=0; j<n; j++){
            sum[i] += rnd.Expon(lambda);
        }
    }
    return sum;
}

double* cauchy_prob_sum(Random rnd, int n, int throws, double mean, double gamma){
    double* sum = new double[throws]{0.};
    for(int i=0; i<throws; i++){
        for(int j=0; j<n; j++){
            sum[i] += rnd.CauchyLorentz(mean, gamma);
        }
    }
    return sum;
}

double* uniform_prob_sum(Random rnd, int n, int throws){
    double* sum = new double[throws]{0.};
    for(int i=0; i<throws; i++){
        for(int j=0; j<n; j++){
            sum[i] += rnd.Rannyu();
        }
    }
    return sum;
}

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int M = 1000000;
   int N[4] = {1,2,10, 100};
   double lambda = 1.;
   double mu = 0.;
   double gamma = 1.;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    //Uniform
    ofstream WriteResults3;
    WriteResults3.open("results_3.dat");
    if (WriteResults3.is_open()){
        double * sum1 = uniform_prob_sum(rnd, N[0], M);
        double * sum2 = uniform_prob_sum(rnd, N[1], M);
        double * sum10 = uniform_prob_sum(rnd, N[2], M);
        double * sum100 = uniform_prob_sum(rnd, N[3], M);
        for(int i=0; i<M; i++) {
            WriteResults3 << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] <<  "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults3.close();

    //Lorentzian
    ofstream WriteResults2;
    WriteResults2.open("results_2.dat");
    if (WriteResults2.is_open()){
        double * sum1 = cauchy_prob_sum(rnd, N[0], M, mu, gamma);
        double * sum2 = cauchy_prob_sum(rnd, N[1], M, mu, gamma);
        double * sum10 = cauchy_prob_sum(rnd, N[2], M, mu, gamma);
        double * sum100 = cauchy_prob_sum(rnd, N[3], M, mu, gamma);
        for(int i=0; i<M; i++) {
            WriteResults2 << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] <<  "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults2.close();

    //Exponential
    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()){
        double * sum1 = expon_prob_sum(rnd, N[0], M, lambda);
        double * sum2 = expon_prob_sum(rnd, N[1], M, lambda);
        double * sum10 = expon_prob_sum(rnd, N[2], M, lambda);
        double * sum100 = expon_prob_sum(rnd, N[3], M, lambda);
        for(int i=0; i<M; i++) {
            WriteResults1 << sum1[i] << " " << sum2[i] << " " << sum10[i] << " " << sum100[i] <<  "\t" << endl;
        }
     } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults1.close();

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
