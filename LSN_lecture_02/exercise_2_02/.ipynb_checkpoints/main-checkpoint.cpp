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

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int M = 10000; //Number of throws
   int n = 3; //Possible directions (both positive and negative)
   int N = 100; //Number of blocks
   int L = M/N; //Number of throws per block;
   double a = 1.; //Lattice step
   double *sum = new double[N]{0};
   int *count = new int[n]{0};
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

   //Discrete case
   for(int j=0; j<N; j++){
       for(int i=0; i<L; i++){
           double r = rnd.Rannyu(-1,1);
           if(r<0){
               count[(int) ((-1)*r * n)] -= a;
           } else if (r>0){
               count[(int) (r * n)] += a;
           }
           sum[i] = sqrt(pow(count[n-3],2)+pow(count[n-2],2)+pow(count[n-1],2));
       }
       count[n] = {0};
   }

   for(int k=0; k<N; k++){
       sum[k] /= N;
       cout << sum[k] << endl;
   }

    ofstream WriteResults1;
    WriteResults1.open("results_1.dat");
    if (WriteResults1.is_open()){
        for(int i=0; i<N; i++){
            WriteResults1 << sum[i] << "\t" << endl;
        }
    } else cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteResults1.close();
   //Continuum case


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
