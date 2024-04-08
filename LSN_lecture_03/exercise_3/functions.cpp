
//
// Created by francesco on 17/03/24.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include "functions.h"
#include "random.h"

using namespace std;

Random initialize(Random rnd, vector<int> seed, int p1, int p2, std::string prime_file, std::string input_file) {
    ifstream Primes(prime_file);
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input(input_file);
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    return rnd;
}





