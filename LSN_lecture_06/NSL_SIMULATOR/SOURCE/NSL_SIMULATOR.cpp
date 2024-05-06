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
#include "system.h"


/*
void equilibration(string method) {
    int nconf = 1;
    int n_part = 50;
    System SYS;
    SYS.initialize_eq(method);
    SYS.initialize_properties_eq(method);

    for (int i = 0; i < SYS.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps(); j++) { //loop over steps in a block
            for (int k = 0; k < n_part; k++){
                SYS.step_eq(k);
                SYS.measure_eq(method);
            }
            if (j % 10 == 0) {
                //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
    }
    return;
}

*/

//La simulazione parte dalla temperatura più elevata perché il sistema, essendo maggiormente disordinato, si trova più vicino
//all'equilibrio, quindi impiega pochi passi per raggiungere l'equilibrio.

//FINIRE IL CODICE DEL MAIN AGGIUNGENDO L'EQUILIBRAZIONE (GUARDANDO I GRAFICI) E STOP

using namespace std;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <value>" << std::endl;
        return 1; // indicate error
    }

    // argv[0] is the name of the program itself
    // argv[1] is the first argument passed by the user

    // Convert the argument to an integer
    string method = argv[1];

    //System SYS;
    //equilibration(method);

    int nconf = 1;
    System SYS;
    SYS.initialize(method);
    SYS.initialize_properties(method);
    SYS.block_reset(0, method);

    for (int i = 0; i < SYS.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps(); j++) { //loop over steps in a block
            SYS.step();
            SYS.measure();
            if (j % 10 == 0) {
                //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        SYS.averages(i + 1, method);
        SYS.block_reset(i + 1, method);
    }
    SYS.finalize(method);

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
