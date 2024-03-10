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
#include <cmath>
#include "random.h"

double error(double * av, double * av2, int n){
    if(n == 0){
        return 0;
    } else {
        return sqrt(av2[n] - pow(av[n],2))/(sqrt(n));
    }
}

double func(double x){
    return (M_PI/2)*cos((M_PI/2)*x);
}

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int M = 10000;
   int N = 100;
   int L = M/N;
   double * mean_integral = new double[L]{0.};
    double * ave = new double[N]{0}; //Arrays for storing the average and the square of the average for each block
    double * ave2 = new double[N]{0}; //Arrays for storing the average and the square of the average for each block
    double * sum_prog = new double[N]{0};
    double * sum2_prog = new double[N]{0};
    double * err_prog = new double[N]{0};
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

    //Block of code for the evaluation of the mean
    for(int i=0; i<N; i++){
        double sum1 = 0.;
        for(int j=0; j<L; j++){
            double sum2 = 0.;
            for(int k=0; k<N; k++){
                double r = rnd.Rannyu();
                sum2 += func(r);
            }
            mean_integral[j] = sum2/N;
            sum1 += mean_integral[j];
            mean_integral[j] = 0;
        }
        ave[i] = double(sum1/L); //Store average values for each block
        ave2[i] = double(pow(ave[i], 2)); //Store square of the average for each block
    }
    for(int k=0; k<N; k++){
        sum_prog[k] = 0.;
        sum2_prog[k] = 0.;
        for(int l=0; l<k+1; l++){
            sum_prog[k] += ave[l];
            sum2_prog[k] += ave2[l];
        }
        sum_prog[k] /= (k+1); //Cumulative average
        sum2_prog[k] /= (k+1); //Cumulative square average
        err_prog[k] = error(sum_prog, sum2_prog, k); //Statistical uncertainty
    }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()){
        for(int i=0; i<N; i++){
            WriteResults1 << sum_prog[i] << " " <<  err_prog[i] << " " << "\t" << endl;
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
