/*
 * CA_GPU.cuh
 *
 *  Created on: 15/apr/2015
 *      Author: knotman
 */

#ifndef CA_GPU_CUH_
#define CA_GPU_CUH_

struct CA_GPU{


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

	//#############################  CA FUNCTIONS (DEVICE)  ##########################
	CA_GPU(){};//default constructor
	void printParameters();




};


//---------------------------------------------------------------------------
void CA_GPU::printParameters()
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
