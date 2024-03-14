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

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int N_experim = 1000000; //Throws
   int N = 100; //Blocks
   int L = N_experim/N; //Throws per block
   double l = 1.; //Length of the needle
   double d = 3.; //Spacing among lines
   int N_hit = 0;
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
        N_hit = 0;
        for(int j=0; j<L; j++){
            double x1 = 0.;
            double x2 = 0.;
            double y1 = 0.;
            double y2 = 0.;
            do {
                y1 = rnd.Rannyu(0,2*d);
                y2 = rnd.Rannyu(0,2*d);
                x1 = rnd.Rannyu(0,2*d);
                x2 = rnd.Rannyu(0,2*d);

            } while(l != sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
            if(((0<x1 && x1<d) && (d<x2 && x2<2*d)) || ((0<x2 && x2<d) && (d<x1 && x1<2*d))){
                N_hit += 1;
            } else if (((0>x1 && x1>(-1)*d) && (0<x2 && x2<d)) || ((0>x2 && x2>(-1)*d) && (0<x1 && x1<d))){
                N_hit += 1;
            }
        }
        cout << N_hit << endl;
        ave[i] = (2*l*L)/(N_hit*d); //Store average values for each block
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
