/*
 * CA_HOST.cuh
 *
 *  Created on: 15/apr/2015
 *      Author: knotman
 */

#ifndef CA_HOST_CUH_
#define CA_HOST_CUH_
#include <iostream>
#include <fstream>
#include <string.h>
#include "utils.cuh"
#include "cuda_error_check.cu"
#include "CA_GPU.cuh"

/**
 * Const char* PATH
 * */
typedef  const char* ccPATH;
/**
 * string PATH
 * */
typedef  std::string sPATH;

struct CA_HOST{
	CA_GPU toBecopied; //HOST handle that contains GPUallocated pointer (needed to call free!)

	sPATH s_data_folder;
	sPATH s_parameters;
	sPATH s_morphology;
	sPATH s_lava_thickness;
	sPATH s_lava_temperature;
	sPATH s_lava_solidified;
	double NODATA_VALUE;
	//#############################  CA PARAMETERS (HOST)  ##########################
	//																				#
	//									 											#
	int h_NR;	//numbers of rows													#
	int h_NC;	//numbers of column													#
	int h_NUMCELLS; //number of cells = h_NR*h_NC									#
	//																				#
	double h_Pclock;	//AC clock [s]												#
	double h_Pc;	//cell side [m]													#
	//																				#
	double h_Pac;	//area of the cell (Pc*Pc) [m^2] 								#
	//																				#
	double h_PTsol;	//temperature of solidification [K]								#
	double h_PTvent;	//temperature of lava at vent [K]							#
	double h_Pr_Tsol; 	//															#
	double h_Pr_Tvent; 	//															#
	//																				#
	double h_a;	// parametro per calcolo Pr - viscosity first parametr				#
	double h_b;	// parametro per calcolo Pr - viscosity second parameter			#
	//																				#
	double h_Phc_Tsol;	//[m]														#
	double h_Phc_Tvent;	//[m]														#
	//																				#
	double h_c;	//parametro per calcolo hc-yield strength first parameter			#
	double h_d;	//parametro per calcolo hc-yield strength seconf parameter			#
	//																				#
	double h_Pcool;	//aderence [m]													#
	double h_Prho;	//density [kg/m^3]												#
	double h_Pepsilon; //emissivity [dimensionless]									#
	double h_Psigma; //Stephen-Boltzamnn constant [J m^-2 s^-1 K^-4]				#
	//																				#
	double h_Pcv; //specific heat [J kg^-1 K^-1]									#
	//																				#
	//###############################################################################



	double* h_sbts; //linearized substates


	//#############################  CA FUNCTIONS (HOST)  ##########################
	CA_HOST();
	void setDataFolderPath(ccPATH dfPath);
	bool loadParameters(ccPATH);

	void simulationInit();
	void evaluatePowerLawParams(double value_sol, double value_vent, double &k1, double &k2);

	bool hostMemoryAllocation();
	void hostMemoryFree();


	void loadSubstateFromFile(ccPATH path, int substate);
	void loadSubstates();

	void printParameters();
	int h_getIndexOfPosition(int x, int y, int substate);
	//#################  CA_GPU CONTROL FUNCTIONS (HOST-SIDE)  ###################
	bool deviceMemoryAllocation(CA_GPU* );	//Device memory allocation
	void deviceMemoryFree(CA_GPU* );		//Device memory free
	CA_GPU* deviceCAGPUInitialization();

	void copyParametersFromCA_HOST_to_CA_GPU(CA_GPU* h_CAGPU);
};


/**
 * It first create an host copy of the GPU structure and copy all the parameter from CA_HOST
 * to this structure. Then the internal buffer (of the host copy of CA_GPU) are allocated ON GPU.
 * The host copy is then copied on an allocated CA:GPU pointer on GPU. All the parameter are copied
 * cleanly and the substates pointer were already pointing to GPU allocated memory!
 * !!!!Should be called ONLY after the CA_HOST initialization is completed!!!!
 * @return a pointer to a fully allocated and initialized GPU structure
 */
CA_GPU* CA_HOST::deviceCAGPUInitialization(){

	copyParametersFromCA_HOST_to_CA_GPU(&toBecopied);
	deviceMemoryAllocation(&toBecopied);
	CUDASAFECALL (cudaMemcpy(toBecopied.d_sbts_current,this->h_sbts,sizeof(double)*h_NUMCELLS,cudaMemcpyHostToDevice));
	//cudaMemcpy(toBecopied.d_sbts_updated,this->h_sbts,sizeof(double)*h_NUMCELLS,cudaMemcpyHostToDevice);//probably unecessary

	CA_GPU* d_pointer;
	CUDASAFECALL(cudaMalloc(&d_pointer,sizeof(CA_GPU)));
	CUDASAFECALL(cudaMemcpy(d_pointer,&toBecopied,sizeof(CA_GPU),cudaMemcpyHostToDevice));

	return d_pointer;
}

bool CA_HOST::deviceMemoryAllocation(CA_GPU* d_CA){
	size_t size = NUMBEROFSUBSTATES*h_NUMCELLS*sizeof(double);
	CUDASAFECALL( cudaMalloc(&d_CA->d_sbts_current,size) );
	CUDASAFECALL( cudaMalloc(&d_CA->d_sbts_updated,size) );
	return true;
}

/**
 * It first free all the GPU array allocated and at
 * the end deallocates the structure itself that lie on GPU
 * @param d_CA
 */
void CA_HOST::deviceMemoryFree(CA_GPU* d_CA){		//Device memory free
	CUDASAFECALL( cudaFree(toBecopied.d_sbts_current) );
	CUDASAFECALL( cudaFree(toBecopied.d_sbts_updated) );
	CUDASAFECALL( cudaFree(d_CA));

}
/**
 * Default Constructor for the HOST-SIDE CA data struct
 * @return void
 */
CA_HOST::CA_HOST(){
	h_NR = 0;
	h_NC = 0;
	h_NUMCELLS=0;
	h_Pclock = 0.0;
	h_Pc = 0.0;
	h_Pac = 0.0;
	h_PTsol = 0.0;
	h_PTvent = 0.0;
	h_Pr_Tsol = 0.0;
	h_Pr_Tvent = 0.0;
	h_a = 0.0;
	h_b = 0.0;
	h_Phc_Tsol = 0.0;
	h_Phc_Tvent = 0.0;
	h_c = 0.0;
	h_d = 0.0;
	h_Pcool = .0;
	h_Prho = 0.0;
	h_Pepsilon = 0.0;
	h_Psigma = 0.0;
	h_Pcv = 0.0;
	h_sbts = NULL;
	//file paths
	s_data_folder="";//root path
	s_parameters="";
	s_lava_thickness="";
	s_lava_temperature="";
	s_lava_solidified="";
	s_morphology="";

	NODATA_VALUE=-9999.0;
}

/**
 * Adjust paths of the input files necessary for the initialization of the CA SCIARAfV2 substates
 * dfPath shoud have a trailing "path separator", '/' in UNIX-like system
 * dfPath should point to a folder that containg files named:
 * 	PARAMETERS.cfg
 * 	MORPHOLOGY.stt
 * 	LAVA_SOLIDIFIED.stt
 * 	LAVA_THICKNESS.stt
 * 	LAVA_TEMPERATURE.stt
 * @param dfPath default the current directory
 */
void CA_HOST::setDataFolderPath(ccPATH dfPath="./"){
	s_data_folder=dfPath;
	s_parameters =s_data_folder+"PARAMETERS.cfg";
	s_morphology = s_data_folder+"MORPHOLOGY.stt";
	s_lava_solidified=s_data_folder+"LAVA_SOLIDIFIED.stt";
	s_lava_thickness=s_data_folder+"LAVA_THICKNESS.stt";
	s_lava_temperature=s_data_folder+"LAVA_TEMPERATURE.stt";
}

/**
 * This function load the parameter from a file specified as parameter.
 * The parameter file should be exactly formatted as below (right column may values may be double numbers
 * except for ncols and nrows that should be integers).
 *
			ncols				448
			nrows				1231
			Pclock				30
			Pc 					10
			PTsol				1143
			PTvent				1360
			Pr_Tsol				0.1
			Pr_Tvent			0.1
			Phc_Tsol			15
			Phc_Tvent			3
			Pcool				3
			Prho				2600
			Pepsilon			0.9
			Psigma				5.68e-8
			Pcv					1150
 *
 * @param path
 * @return a bool value which value tell if parsing succeded.
 */
bool CA_HOST::loadParameters(ccPATH path){
	FILE *file;
	char str[255];
	const char ncols_str[] 		= "ncols";
	const char nrows_str[] 		= "nrows";
	const char Pclock_str[] 	= "Pclock";
	const char Pc_str[] 		= "Pc";
	const char PTsol_str[] 		= "PTsol";
	const char PTvent_str[] 	= "PTvent";
	const char Pr_Tsol_str[] 	= "Pr_Tsol";
	const char Pr_Tvent_str[] 	= "Pr_Tvent";
	const char Phc_Tsol_str[] 	= "Phc_Tsol";
	const char Phc_Tvent_str[] 	= "Phc_Tvent";
	const char Pcool_str[] 		= "Pcool";
	const char Prho_str[] 		= "Prho";
	const char Pepsilon_str[] 	= "Pepsilon";
	const char Psigma_str[] 	= "Psigma";
	const char Pcv_str[] 		= "Pcv";

	if (( file = fopen(path, "r"))==NULL)
	{
		fprintf(stderr,"Cannot open file parameters.\nEXITING");
		exit(1);
	}

	//ncols
	fscanf(file,"%s",&str);
	if (strcmp(str, ncols_str)){
		fprintf(stderr,"Error ncols.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_NC = atoi(str);

	//nrows
	fscanf(file,"%s",&str);
	if (strcmp(str, nrows_str)){
		fprintf(stderr,"Error nrows.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_NR = atoi(str);

	//Pclock
	fscanf(file,"%s",&str);
	if (strcmp(str, Pclock_str)){
		fprintf(stderr,"Error Pclock.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pclock = atof(str);

	//Pc = cell_size
	fscanf(file,"%s",&str);
	if (strcmp(str, Pc_str)){
		fprintf(stderr,"Error Pc.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pc = atof(str);

	h_Pac=h_Pc*h_Pc;

	//PTsol
	fscanf(file,"%s",&str);
	if (strcmp(str, PTsol_str)){
		fprintf(stderr,"Error PTsol.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_PTsol = atof(str);

	//PTvent
	fscanf(file,"%s",&str);
	if (strcmp(str, PTvent_str)){
		fprintf(stderr,"Error PTvent.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_PTvent = atof(str);

	//Pr_Tsol
	fscanf(file,"%s",&str);
	if (strcmp(str, Pr_Tsol_str)){
		fprintf(stderr,"Error Pr_Tsol.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pr_Tsol = atof(str);

	//Pr_Tvent
	fscanf(file,"%s",&str);
	if (strcmp(str, Pr_Tvent_str)){
		fprintf(stderr,"Error Pr_Tvent.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pr_Tvent = atof(str);

	//Phc_Tsol
	fscanf(file,"%s",&str);
	if (strcmp(str, Phc_Tsol_str)){
		fprintf(stderr,"Error Phc_Tsol.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Phc_Tsol = atof(str);

	//Phc_Tvent
	fscanf(file,"%s",&str);
	if (strcmp(str, Phc_Tvent_str)){
		fprintf(stderr,"Error Phc_Tvent.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Phc_Tvent = atof(str);

	//Pcool
	fscanf(file,"%s",&str);
	if (strcmp(str, Pcool_str)){
		fprintf(stderr,"Error Pcool.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pcool = atof(str);

	//Prho
	fscanf(file,"%s",&str);
	if (strcmp(str, Prho_str)){
		fprintf(stderr,"Error Prho.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Prho = atof(str);

	//Pepsilon
	fscanf(file,"%s",&str);
	if (strcmp(str, Pepsilon_str)){
		fprintf(stderr,"Error Pepsilon.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pepsilon = atof(str);

	//Psigma
	fscanf(file,"%s",&str);
	if (strcmp(str, Psigma_str)){
		fprintf(stderr,"Error Psigma.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Psigma = atof(str);

	//Pcv
	fscanf(file,"%s",&str);
	if (strcmp(str, Pcv_str)){
		fprintf(stderr,"Error Pcv.\n");
		return false;
	}
	fscanf(file,"%s",&str);
	h_Pcv = atof(str);

	fclose(file);

	return true;

}//load parameters

void CA_HOST::evaluatePowerLawParams(double value_sol, double value_vent, double &k1, double &k2){

	k2 = ( log10(value_vent) - log10(value_sol) ) / (h_PTvent - h_PTsol) ;
	k1 = log10(value_sol) - k2*(h_PTsol);
}

/**
 * Performs all the preliminary operation required
 * to finalize the configuration of the simulation
 * -compute the a,b,c,d parameters
 * - allocate memory for substates
 */
void CA_HOST::simulationInit(){
	//compute the a,b,c,d parameters
	evaluatePowerLawParams(h_Pr_Tsol, h_Pr_Tvent, h_a, h_b);
	evaluatePowerLawParams(h_Phc_Tsol, h_Phc_Tvent, h_c, h_d);
	h_NUMCELLS=h_NC*h_NR;

	//allocate the memory for the substates and other CA structures

	bool go= hostMemoryAllocation();
	if(!go){
		fatalErrorExit("ALLOCATION ERROR");
	}
}

/**
 * Allocate
 * @return Was allocation succesfull? Did calloc return a valid heap pointer?
 */
bool CA_HOST::hostMemoryAllocation(){
	//linearized substates.
	h_sbts = (double *) calloc(NUMBEROFSUBSTATES*h_NR*h_NC, sizeof(double));
	if(h_sbts)
		return true;
	//if something goes wrong with the memory allocation
	return false;
}

/**
 * Release memory used for the substates of the CA host side
 * Set the pointer to NULL
 */
void CA_HOST::hostMemoryFree(){
	free(h_sbts);
	h_sbts=NULL;
}


/**
 * Load the substate, reading the file from the path parameter
 * @param path
 * @param substate
 */
void CA_HOST::loadSubstateFromFile(ccPATH path, int substate){

	std::ifstream InFile(path);

	if (!InFile){
		InFile.close();
		std::string error = "Error opening the file ";
		error+= path;
		fatalErrorExit(error.c_str());
		//implicit exit here calling fatalErrorExit; Execution stops here!
	}

	double pointValue;
	for(int row=0; row<h_NR;row++) {
		for(int col=0; col<h_NC;col++){

			InFile >> pointValue; /* reads a single value for a cell */

			if (pointValue == NODATA_VALUE)
				pointValue = 0;

			h_sbts[h_getIndexOfPosition(row, col, substate)] = pointValue;
		}
	}
	InFile.close();
}

/**
 * Read and initialize the substates from file
 * There is no need to initialize all the substates
 * (since starting from scratch a simulation, there
 * is no lava in the cellular space and hence no
 * need for thinkness,  * temperature and solidified
 * lava, and any flows).  It really make sense if want
 *  to implement the hot start of a simulation
 *  (from a previous incomplete simulation).
 *  One can think to save the partial result in file
 *  for substates temperature thickness and solidified
 *  lava and restore the simulation.
 */
void CA_HOST::loadSubstates(){
	loadSubstateFromFile(s_morphology.c_str(),QUOTE);
}

//---------------------------------------------------------------------------
void CA_HOST::printParameters()
{
	printf("---------------------------------------------\n");
	printf("Paramater		Value\n");
	printf("---------------------------------------------\n");
	printf("ncols			%u\n",  h_NC);
	printf("nrows			%u\n",  h_NR);
	printf("Pclock			%f\n",  h_Pclock);
	printf("PTsol			%f\n",  h_PTsol);
	printf("PTvent			%f\n",  h_PTvent);
	printf("Pr(Tsol)		%f\n",  h_Pr_Tsol);
	printf("Pr(Tvent)		%f\n",  h_Pr_Tvent);
	printf("Phc(Tsol)		%f\n",  h_Phc_Tsol);
	printf("Phc(Tvent)		%f\n",  h_Phc_Tvent);
	printf("Pcool			%f\n",  h_Pcool);
	printf("Prho			%f\n",  h_Prho);
	printf("Pepsilon		%f\n",  h_Pepsilon);
	printf("Psigma			%e\n",  h_Psigma);
	printf("Pcv			%f\n",  h_Pcv);
	printf("a			%f\n",  h_a);
	printf("b			%f\n",  h_b);
	printf("c			%f\n",  h_c);
	printf("d			%f\n",  h_d);
	printf("---------------------------------------------\n");
}

/**
 * The linearized index computation works intuitively as follows:
 * Fist compleate substate 2d matrices have to be jumped (substate*NUMCELLS = substate*(h_NC*h_NR) )
 * Then we can retrieve the index inside the h_NC*h_NR cells of the right substate matrix
 * @param x
 * @param y
 * @param substate
 * @return the 2D row-major linearized index of the corrensponding 2D array location (x,y)
 */
__host__  int CA_HOST::h_getIndexOfPosition(int x, int y, int substate){
	return ( (h_NUMCELLS * substate)  +   (x * h_NC) + y );
}

void CA_HOST::copyParametersFromCA_HOST_to_CA_GPU(CA_GPU* h_CAGPU){
	if(!h_CAGPU)
		fatalErrorExit("NULL pointer to CA_GPU host structure");

	h_CAGPU->d_NR 				= this->h_NR;
	h_CAGPU->d_NC 				= this->h_NC;
	h_CAGPU->d_NUMCELLS 		= this->h_NUMCELLS;
	h_CAGPU->d_Pclock 			= this->h_Pclock;
	h_CAGPU->d_Pc 				= this->h_Pc;
	h_CAGPU->d_Pac 				= this->h_Pac;
	h_CAGPU->d_PTsol 			= this->h_PTsol;
	h_CAGPU->d_PTvent  			= this->h_PTvent;
	h_CAGPU->d_Pr_Tsol 			= this->h_Pr_Tsol;
	h_CAGPU->d_Pr_Tvent 		= this->h_Pr_Tvent;
	h_CAGPU->d_a 				= this->h_a;
	h_CAGPU->d_b 				= this->h_b;
	h_CAGPU->d_Phc_Tsol 		= this->h_Phc_Tsol;
	h_CAGPU->d_Phc_Tvent 		= this->h_Phc_Tvent;
	h_CAGPU->d_c 				= this->h_c;
	h_CAGPU->d_d 				= this->h_d;
	h_CAGPU->d_Pcool 			= this->h_Pcool;
	h_CAGPU->d_Prho 			= this->h_Prho;
	h_CAGPU->d_Pepsilon 		= this->h_Pepsilon;
	h_CAGPU->d_Psigma 			= this->h_Psigma;
	h_CAGPU->d_Pcv 				= this->h_Pcv;
}

#endif /* CA_HOST_CUH_ */
