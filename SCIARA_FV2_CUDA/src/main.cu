#include <stdio.h>
#include <stdlib.h>

//####  OTHER INCLUDES  #########
#include "CA_HOST.cuh"//		#
#include "utils.cuh"

//###############################

//####  GLOBAL VARIABLES  #########
CA_HOST h_CA;
CA_GPU* d_CA;

dim3 dimBlock;
dim3 dimGrid;
//###############################
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
	//d_CA->printSubstate(0);
	printf("NUM of vents: %u ",d_CA->numVents);
	for(int i =0;i<d_CA->numVents;i++){
		printf("(%u %u)\n",d_CA->coordVents[i*2],d_CA->coordVents[i*2+1]);
	}

	for(int i =0;i<d_CA->numVents;i++){
		for(int j=0;j<d_CA->emissionRate_size;j++){
			printf("\tEMISSION RATE vent %u time %u -> %f\n",i,j,d_CA->emissionRates[d_CA->emissionRate_size*i+j]);
		}
	}
}

/**
 * Lava emission from vents.
 * This kernel should be launched with wust one block in 1D
 * (on X dimension) and with the number of threads equals to the number of
 * vents (parameter numVents)
 * @param d_CA
 */
__global__ void emitLavaFromVents(CA_GPU* d_CA){

	if(blockIdx.x==0 && blockIdx.y==0){
		if(threadIdx.x < d_CA->numVents && threadIdx.y==0){
			d_CA->emitLavaFromVent(threadIdx.x);
		}
	}
}

/**
 * Temperature initialization
 * @param d_CA
 */
__global__ void temperatureInitialization(CA_GPU* d_CA){
	d_CA->cellTemperatureInitialize();
}

/**
 * Flows Computation
 * @param d_CA
 */
__global__ void computeFlows(CA_GPU* d_CA){
	d_CA->empiricalFlows();
}

/**
 * Flows reduction
 * @param d_CA
 */
__global__ void reduceFlows(CA_GPU* d_CA){
	d_CA->distribuiteFlows();
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
unsigned int nVents;
void globalTransitionFunction(){
	////kernel launch parameters settings
	dimBlock.x=8;
	dimBlock.y=8;
	computeKernelLaunchParameter(dimBlock.x,dimBlock.y,h_CA.getNr(),h_CA.getNc(),dimGrid);


	//for(int i=0;i<100;i++){
	emitLavaFromVents<<<1,nVents>>>(d_CA);
	temperatureInitialization<<<dimGrid,dimBlock>>>(d_CA);
	computeFlows<<<dimGrid,dimBlock>>>(d_CA);
	//}

}



int main ( int argc, char *argv[] ){
	cudaDeviceReset();
	hostInitialize(argc,argv);
	//configure CA HOST
	h_CA.simulationInit();
	nVents = h_CA.getNumVents();

	h_CA.loadSubstates();
	h_CA.printParameters();

	d_CA=h_CA.deviceCAGPUInitialization();




	/*
	 * GLOBAL TRANSITION FUNCTION ON GPU
	 */

	globalTransitionFunction();





	h_CA.copyBackFromGPU(d_CA);
	h_CA.saveSubstatesOnFile(h_CA.getDataFolder()+"/output/");

	//host initialization and configuration completed
	h_CA.deviceMemoryFree(d_CA);
	//free CA_HOST memory
	h_CA.hostMemoryFree();
	printf("SIMULATION ENDED\n");
}
