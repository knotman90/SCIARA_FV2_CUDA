#include <stdio.h>
#include <stdlib.h>

//####  OTHER INCLUDES  #########
#include "CA_HOST.cuh"//		#
#include "utils.cuh"
#include "timeCounter.h"

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
	if(!h_CA.loadParameters()){
		fatalErrorExit("FAILED LOADING CONFIGURATION PARAMETERS");
	}
}



//##### TRANSITION FUNCTION KERNELS ######

__global__ void printSubstateG(CA_GPU* d_CA, int substate){
	d_CA->printSubstate(substate);
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
 * It computes the new values of thicknes based on the previous step of
 * computation, hence the calculation of the outflows
 * It also compute the new lava temperature and handle the lava solidification
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

/** NOT A KERNEL!
 * Copy matrix.
 * Hard Swap the content  current and updated matrix
 * Lanches a kernel in which thread explicitly copy the matrices from
 * the current to the updated version
 * @param d_CA
 */
// void copyMatricesMemCpyDevToDev(CA_GPU* d_CA){
//	h_CA->copyMatricesMemCpyDevToDev();
//}

/**
 * Copy matrix.
 * Hard Swap the content  current and updated matrix
 * Lanches a kernel in which thread explicitly copy the matrices from
 * the current to the updated version
 * @param d_CA
 */
__global__ void copyMatrices(CA_GPU* d_CA){
	d_CA->swapMatrices();
}

__global__ void updateMinRect(CA_GPU* d_CA){
	d_CA->copyNewToCurrentAdaptiveGrid();
}

void printfAdaptoveGrid(){
	for(int i=0;i<ADAPTIVEGRID_SIZE;i++){
		printf("%i ",h_CA.h_d_adaptive_grid[i]);
	}
	printf("\n");
}

//#######################################
unsigned int nVents;
void globalTransitionFunction(){
	/*To print adaptive grid
	 * cudaDeviceSynchronize();
		printf("START ");
		printfAdaptoveGrid();//delete just for debug
		cudaDeviceSynchronize();
	 *
	 */
	////kernel launch parameters settings
	int COLS,ROWS;
	dimBlock.x=8;
	dimBlock.y=8;

	computeKernelLaunchParameter_plus(dimBlock,h_CA.h_d_adaptive_grid,1,dimGrid);

#pragma unroll
	for(int i=0;i<h_CA.getNSteps();i++){
		//lava emission from vents
		emitLavaFromVents<<<1,nVents>>>(d_CA);
		//transformation of the temperature in temp*thickness
		temperatureInitialization<<<dimGrid,dimBlock>>>(d_CA);
		//outflows computations. Each cells compute flow out to
		computeFlows<<<dimGrid,dimBlock>>>(d_CA);
		computeKernelLaunchParameter_plus(dimBlock,h_CA.h_d_adaptive_grid,3,dimGrid);//fat adaptivegrid to handle expansion
		reduceFlows<<<dimGrid,dimBlock>>>(d_CA);
		updateMinRect<<<1,4>>>(d_CA);
		computeKernelLaunchParameter_plus(dimBlock,h_CA.h_d_adaptive_grid,1,dimGrid);
		//h_CA.copyMatricesMemCpyDevToDev();
		copyMatrices<<<dimGrid,dimBlock>>>(d_CA);
	}

}



__global__ void globalTransitionFunctionGPU(uint NR,uint NC, uint nSteps,uint nVents,CA_GPU* d_CA){
	dim3 dimBlock;
	dim3 dimGrid;
	dimBlock.x=8;
	dimBlock.y=8;
	computeKernelLaunchParameter(dimBlock.x,dimBlock.y,NR,NC,dimGrid);
#pragma unroll
	for(int s=0;s<nSteps;s++){
		emitLavaFromVents<<<1,nVents>>>(d_CA);
		temperatureInitialization<<<dimGrid,dimBlock>>>(d_CA);
		computeFlows<<<dimGrid,dimBlock>>>(d_CA);
		reduceFlows<<<dimGrid,dimBlock>>>(d_CA);
		copyMatrices<<<dimGrid,dimBlock>>>(d_CA);
	}
}

void gtf_kernel(){
	globalTransitionFunctionGPU<<<1,1>>>(h_CA.getNr(),h_CA.getNc(),h_CA.getNSteps(), h_CA.getNumVents(),d_CA);
	cudaDeviceSynchronize();
}

__global__ void testUnified(CA_GPU* d_CA){
	printf("%i\n",d_CA->h_d_adaptive_grid[0]);
}

int main ( int argc, char *argv[] ){
	cudaDeviceReset();

	hostInitialize(argc,argv);
	//configure CA HOST
	h_CA.simulationInit();
	nVents = h_CA.getNumVents();

	h_CA.loadSubstates();
	//	h_CA.printParameters();

	d_CA=h_CA.deviceCAGPUInitialization();




	/*
	 * GLOBAL TRANSITION FUNCTION ON GPU
	 */
	//globalTransitionFunctionGPU<<<1,1>>>(h_CA.getNr(),h_CA.getNc(),h_CA.getNSteps(), h_CA.getNumVents(),d_CA);
	//gtf_kernel();
	auto elapsedTime= tim::measure<>::execution(globalTransitionFunction); //CPU managed
	//auto elapsedTime= tim::measure<>::execution(gtf_kernel); //GPU managed
	//globalTransitionFunction();





	h_CA.copyBackFromGPU(d_CA);
	h_CA.saveSubstatesOnFile(h_CA.getDataFolder()+"/output/");

	//host initialization and configuration completed
	h_CA.deviceMemoryFree(d_CA);
	//free CA_HOST memory
	h_CA.hostMemoryFree();
	//printf("SIMULATION ENDED in %i ms\n",elapsedTime);
	printf("SIMULATION ENDED\n");
}
