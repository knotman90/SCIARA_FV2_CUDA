#include <stdio.h>
#include <stdlib.h>

//####  OTHER INCLUDES  #########
#include "CA_HOST.h"//		#
#include "utils.h"
//###############################
CA_HOST h_CA;
/**
 * Parse the command line argument and consequently sets
 * the appropriate file paths for the parameters of the simulation.
 * @param argc
 * @param argv
 */
void hostInitialize(int argc, char *argv[]){
	//get options and file paths from command arguments
	struct CommandLine cmd;
	cmd.parseArgv(argc,argv);
	//configure hosts simulation
	h_CA.setDataFolderPath(cmd._load_path);
	h_CA.loadParameters(h_CA.s_parameters.c_str());
}


int main ( int argc, char *argv[] ){
	hostInitialize(argc,argv);
	//configure CA HOST
	h_CA.simulationInit();
	//allocate CA_HOST memory
	h_CA.hostMemoryAllocation();
	h_CA.loadSubstates();
	h_CA.printParameters();

	//free CA_HOST memory
	h_CA.hostMemoryFree();
	printf("SIMULATION ENDED\n");
}
