#include <stdio.h>
#include <stdlib.h>

//####  OTHER INCLUDES  #########
#include "CA_HOST.cuh"//		#
#include "utils.cuh"

//###############################
CA_HOST h_CA;
CA_GPU* d_CA;
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
	h_CA.loadParameters();
}



//##### TRANSITION FUNCTION KERNELS ######

__global__ void printSubstateG(CA_GPU* d_CA, int substate){
	d_CA->printSubstate(0);
}

/**
 * Lava emission from vents
 * @param d_CA
 */
__global__ void emitLavaFromVents(CA_GPU* d_CA){

}

/**
 * Temperature initialization
 * @param d_CA
 */
__global__ void temperatureInitialization(CA_GPU* d_CA){

}

/**
 * Flows Computation
 * @param d_CA
 */
__global__ void computeFlows(CA_GPU* d_CA){

}

/**
 * Flows reduction
 * @param d_CA
 */
__global__ void reduceFlows(CA_GPU* d_CA){

}

/**
 * Temperture Update reduction
 * @param d_CA
 */
__global__ void temperatureUpdate(CA_GPU* d_CA){

}

/**
 * Copy matrix. Hard Swap the content  current and updated matrix
 * @param d_CA
 */
__global__ void copyMatrices(CA_GPU* d_CA){

}

//#######################################
void globalTransitionFunction(){

}



int main ( int argc, char *argv[] ){
	cudaDeviceReset();
	hostInitialize(argc,argv);
	//configure CA HOST
	h_CA.simulationInit();

	h_CA.loadSubstates();
	h_CA.printParameters();

	d_CA=h_CA.deviceCAGPUInitialization();

/*
 * GLOBAL TRANSITION FUNCTION ON GPU
 */

	h_CA.copyBackFromGPU(d_CA);
	h_CA.saveSubstatesOnFile(h_CA.getDataFolder()+"/output/");

	//host initialization and configuration completed
	h_CA.deviceMemoryFree(d_CA);
	//free CA_HOST memory
	h_CA.hostMemoryFree();
	printf("SIMULATION ENDED\n");
}
