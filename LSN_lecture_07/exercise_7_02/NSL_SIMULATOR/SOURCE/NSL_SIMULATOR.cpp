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

using namespace std;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <value>" << std::endl;
        return 1; // indicate error
    }

    // argv[0] is the name of the program itself
    // argv[1] is the first argument passed by the user

    // Convert the argument to an integer
    string phase = argv[1];


    //Equilibration phase
    int eq_steps = 1000;
    int nconf = 1;
    System SYS;
    SYS.initialize(phase);
    SYS.initialize_properties(phase);
    SYS.block_reset(0, phase);

    for (int i = 0; i < eq_steps; i++) { //equilibration steps
        SYS.step(phase);
    }
    SYS.write_configuration(phase);

    //RESTART THE SIMULATION
    SYS.block_reset(0, phase);
    SYS.read_configuration(phase);
    SYS.initialize_velocities(phase);
    for (int i = 0; i < SYS.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps(); j++) { //loop over steps in a block
            SYS.step(phase);
            SYS.measure(phase);
            if (j % 10 == 0) {
//              SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        //SYS.averages(i + 1, phase);
        SYS.block_reset(i + 1, phase);
    }

    SYS.finalize(phase);

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
