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
        std::cerr << "Usage: " << argv[0] << " <value>" << argv[1] << "<value>" << std::endl;
        return 1; // indicate error
    }

    // argv[0] is the name of the program itself
    // argv[1] is the first argument passed by the user
    // argv[2] is the second argument passed by the user


    // Convert the argument to an integer
    string phase = argv[1];
    string sim_type = argv[2];

    bool breaking = false;
    int nconf = 1;
    System SYS;
    SYS.initialize(phase, sim_type);
    SYS.initialize_properties(phase, sim_type);
    SYS.block_reset(0, phase, sim_type);

    for (int i = 0; i < SYS.get_nbl() && !breaking; i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps() && !breaking; j++) { //loop over steps in a block
            if (!SYS.get_restart()) {
                breaking = true;
                SYS.write_configuration(phase, sim_type);
            } else {
                SYS.step(phase);
                SYS.measure_temp();
            }
        }
    }

    //AGGIUNGERE PRIMI STEP DI EQUILIBRAZIONE COME IN 7.2 e 6 !!!!!
    //SISTEMARE PRESSIONE, VIENE LEGGERMENTE DIVERSA DA CORTI
    //RESTART THE SIMULATION
    SYS.block_reset(0, phase, sim_type);
    SYS.read_configuration(phase, sim_type);
    SYS.initialize_velocities(phase, sim_type);
    for (int i = 0; i < SYS.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps(); j++) { //loop over steps in a block
            SYS.step_restart(phase);
            SYS.measure();
            if (j % 10 == 0) {
//              SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        SYS.averages(i + 1, phase, sim_type);
        SYS.block_reset(i + 1, phase, sim_type);
    }

    SYS.finalize(phase, sim_type);

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
