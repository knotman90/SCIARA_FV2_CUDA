/*
 * CA_GPU.cuh
 *
 *  Created on: 15/apr/2015
 *      Author: knotman
 */

#ifndef CA_GPU_CUH_
#define CA_GPU_CUH_
class CA_GPU{
	friend class CA_HOST;
private:
	//#############################  CA VARIABLES (DEVICE)  ##########################
	double d_sim_elapsed_time; //tempo trascorso dall'inizio della simulazione [s]

	//#############################  CA PARAMETERS (DEVICE)  ##########################
	//																				#
	//									 											#
	int d_NR;	//numbers of rows													#
	int d_NC;	//numbers of column													#
	int d_NUMCELLS; //number of cells = d_NR*d_NC									#
	//																				#
	double d_Pclock;	//AC clock [s]												#
	double d_Pc;	//cell side [m]													#
	//																				#
	double d_Pac;	//area of the cell (Pc*Pc) [m^2] 								#
	//																				#
	double d_PTsol;	//temperature of solidification [K]								#
	double d_PTvent;	//temperature of lava at vent [K]							#
	double d_Pr_Tsol; 	//															#
	double d_Pr_Tvent; 	//															#
	//																				#
	double d_a;	// parametro per calcolo Pr - viscosity first parametr				#
	double d_b;	// parametro per calcolo Pr - viscosity second parameter			#
	//																				#
	double d_Phc_Tsol;	//[m]														#
	double d_Phc_Tvent;	//[m]														#
	//																				#
	double d_c;	//parametro per calcolo hc-yield strength first parameter			#
	double d_d;	//parametro per calcolo hc-yield strength seconf parameter			#
	//																				#
	double d_Pcool;	//aderence [m]													#
	double d_Prho;	//density [kg/m^3]												#
	double d_Pepsilon; //emissivity [dimensionless]									#
	double d_Psigma; //Stephen-Boltzamnn constant [J m^-2 s^-1 K^-4]				#
	//																				#
	double d_Pcv; //specific heat [J kg^-1 K^-1]									#
	//																				#
	//###############################################################################

	//double matrix strategy for CA substated evolution
	double* d_sbts_current; //linearized substates
	double* d_sbts_updated; //linearized substates

public:
	//### VENTS MANAGEMENT AND EMISSION RATES ###########################
	unsigned int emission_time;
	unsigned int numVents;
	unsigned int emissionRate_size;
	int *coordVents;//dim=2*numVents
	double * emissionRates;//dim=numVents*sizeEmissionRate;
	//###################################################################

	//#############################  CA FUNCTIONS (DEVICE)  ##########################
	CA_GPU(){};//default constructor
	__host__ void printParameters();


	//### KERNELS AND GPU CODE FOR THE TRANSITION FUNCTION #####
	__device__ int d_getIdx(int x, int y, int substate)	{
		return ( (d_NUMCELLS * substate)  +   (x * d_NC) + y );
	}

	__device__ void printSubstate(int substate);

	__device__ double ventThickness(unsigned int vent);
	__device__ void emitLavaFromVent(unsigned int vent);

};

__device__ double CA_GPU::ventThickness(unsigned int vent){
	unsigned int i;
	i = (unsigned int) (d_sim_elapsed_time / emission_time);
	if (i >= emissionRate_size)
		return 0;
	else
		return emissionRates[vent*emissionRate_size+i] / d_Pac * d_Pclock;
}


__device__ void CA_GPU::emitLavaFromVent(unsigned int vent){
	double emitted_lava = 0;
	//slt = lava thickness current

	//emit lava
	emitted_lava = ventThickness(vent);
	if (emitted_lava > 0) {
		int x = get_X_LinearIdxVentCoord(vent);
		int y = get_Y_LinearIdxVentCoord(vent);
		int cellIdx = d_getIdx(x,y,THICKNESS);
		d_sbts_current[cellIdx] += emitted_lava;
		d_sbts_updated[cellIdx] = d_sbts_current[cellIdx];

		d_sbts_current[d_getIdx(x,y,TEMPERATURE)] = d_PTvent;
		d_sbts_updated[d_getIdx(x,y,TEMPERATURE)] = d_PTvent;
	}

}


__device__ void CA_GPU::printSubstate(int substate){
	for(int i=0;i<10;i++){
		for (int j = 0; j < 10; j++) {
			printf("%.3f ",d_sbts_current[d_getIdx(i,j,0)]);
		}
		printf("\n");
	}
}
//---------------------------------------------------------------------------
__host__ void CA_GPU::printParameters()
{
	printf("---------------------------------------------\n");
	printf("Paramater		Value\n");
	printf("---------------------------------------------\n");
	printf("Pclock			%f\n",  d_Pclock);
	printf("PTsol			%f\n",  d_PTsol);
	printf("PTvent			%f\n",  d_PTvent);
	printf("Pr(Tsol)		%f\n",  d_Pr_Tsol);
	printf("Pr(Tvent)		%f\n",  d_Pr_Tvent);
	printf("Phc(Tsol)		%f\n",  d_Phc_Tsol);
	printf("Phc(Tvent)		%f\n",  d_Phc_Tvent);
	printf("Pcool			%f\n",  d_Pcool);
	printf("Prho			%f\n",  d_Prho);
	printf("Pepsilon		%f\n",  d_Pepsilon);
	printf("Psigma			%e\n",  d_Psigma);
	printf("Pcv			%f\n",  d_Pcv);
	printf("---------------------------------------------\n");
}

#endif /* CA_GPU_CUH_ */
