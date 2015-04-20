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
	__inline__
	__device__ int getMooreNeighIdx_X(int row, int neighIdx);

	__inline__
	__device__ int getMooreNeighIdx_Y(int col, int neighIdx);

	__inline__
	__device__ int d_getIdx(int x, int y, int substate)	{
		return ( (d_NUMCELLS * substate)  +   (x * d_NC) + y );
	}

	__inline
	__device__ bool isWithinCABounds(int row, int col);

	__inline__
	__device__ void printSubstate(int substate);

	//vent emission--------------------------------------------
	__inline__
	__device__ double ventThickness(unsigned int vent);

	__inline__
	__device__ void emitLavaFromVent(unsigned int vent);
	//----------------------------------------------------------

	__inline__
	__device__ void cellTemperatureInitialize();


	__device__ void empiricalFlows();

	__device__ double __inline__ PowerLaw(double k1, double k2, double T){
		return exp10(k1 + k2*T);
	}
	__device__ void distribuiteFlows();



};//end definition of CA_GPU

__device__ void CA_GPU::distribuiteFlows(){
	//get cell coordinates
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if(isWithinCABounds(row,col)){
		double new_tick = d_sbts_updated[d_getIdx(row,col,THICKNESS)];
		//flown(x,y) has to be added comes from the sud cell
		//flows(neighbor1) has to be substracted because is the amount of lava that the current cell
		//transfer to the upper cell
//same hold for the other pair of indices


		//NORD
		new_tick-= d_sbts_updated[d_getIdx(row,col,FLOWN)] -
				d_sbts_updated[d_getIdx(getMooreNeighIdx_X(row,1),getMooreNeighIdx_Y(col,1),FLOWS)];
	}
}

/**
 * Return the X coordinate of the neigh number neighIdx
         5 | 1 | 8
        ---|---|---
         2 | 0 | 3
        ---|---|---
         6 | 4 | 7
 */
__inline__
__device__ int CA_GPU::getMooreNeighIdx_X(int row, int neighIdx){
	if(neighIdx==0 || neighIdx==2 || neighIdx == 3)
		return row;

	if(neighIdx==1 || neighIdx==5 || neighIdx == 8)
		return row-1;

	//implicit else if(neighIdx==4 || neighIdx==6 || neighIdx == 7)
	return row+1;
}


/**
 * Return the Y coordinate of the neigh number neighIdx
         5 | 1 | 8
        ---|---|---
         2 | 0 | 3
        ---|---|---
         6 | 4 | 7
 */
__inline__
__device__ int CA_GPU::getMooreNeighIdx_Y(int col, int neighIdx){
	if(neighIdx==0 || neighIdx==1 || neighIdx == 4)
		return col;

	if(neighIdx==2 || neighIdx==5 || neighIdx == 6)
		return col-1;

	//implicit else if(neighIdx==3 || neighIdx==7 || neighIdx == 8)
	return col+1;

}




__inline
__device__ bool CA_GPU::isWithinCABounds(int row, int col){

	return
			(row>=0 && row <= (d_NR-1)) &&
			col>=0 && col <= (d_NC-1);
}


__inline__
__device__ double CA_GPU::ventThickness(unsigned int vent){
	unsigned int i;
	i = (unsigned int) (d_sim_elapsed_time / emission_time);
	if (i >= emissionRate_size)
		return 0;
	else
		return emissionRates[vent*emissionRate_size+i] / d_Pac * d_Pclock;
}


__inline__
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


__inline__
__device__ void CA_GPU::cellTemperatureInitialize(){
	//get cell coordinates
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if(isWithinCABounds(row,col)){

		//current cell temperature
		double temp = d_sbts_current[d_getIdx(row,col,TEMPERATURE)];

		if(temp>=0)
			d_sbts_updated[d_getIdx(row,col,TEMPERATURE)] = temp*temp;

		return;
	}

}

__device__ void CA_GPU::empiricalFlows(){
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	if(isWithinCABounds(row,col)){

		if (d_sbts_current[d_getIdx(row,col,THICKNESS)] > 0) {

			bool n_eliminated[MOORE_NEIGHBORS];
			double z[MOORE_NEIGHBORS];
			double h[MOORE_NEIGHBORS];
			double H[MOORE_NEIGHBORS];
			double theta[MOORE_NEIGHBORS];
			double w[MOORE_NEIGHBORS];		//Distances between central and adjecent cells
			double Pr[MOORE_NEIGHBORS];		//Relaxation rate array
			bool loop;
			int counter,i,j;
			double avg;
			double _w,_Pr,hc;
			_w 	= d_Pc;
			_Pr = PowerLaw(d_a, d_b, d_sbts_current[d_getIdx(row,col,TEMPERATURE)]);
			hc 	= PowerLaw(d_c, d_d, d_sbts_current[d_getIdx(row,col,TEMPERATURE)]);

#pragma loop unroll
			for ( i = 0; i < MOORE_NEIGHBORS; i++) {

				h[i] = d_sbts_current[d_getIdx(getMooreNeighIdx_X(row,i),getMooreNeighIdx_Y(col,i),THICKNESS)] ;
				H[i]  = theta[i] = 0;
				w[i] = _w;
				Pr[i] = _Pr;

				if (i < VON_NEUMANN_NEIGHBORS)
					z[i] = d_sbts_current[d_getIdx(getMooreNeighIdx_X(row,i),getMooreNeighIdx_Y(col,i),ALTITUDE)] ;
				else
					//val - (val-valNeigh_I)/rad2
					z[i]= d_sbts_current[d_getIdx(row,col,ALTITUDE)] -
					(
							d_sbts_current[d_getIdx(row,col,ALTITUDE)] -
							d_sbts_current[d_getIdx(getMooreNeighIdx_X(row,i),getMooreNeighIdx_Y(col,i),ALTITUDE)]

					)/RAD2;
			}//for

			H[0] = z[0];
			n_eliminated[0]=true;
#pragma loop unroll
			for ( i = 1; i < MOORE_NEIGHBORS; i++){
				/*//Donato's code is commented here while Giuseppe's not
				//Stick with official Donato version
				if ((z[0] > z[i]) && (h[0] >= h[i]))
				{
					H[i] = z[i];
					theta[i] = atan( (z[0] - z[i]) / w[i] );
					n_eliminated[i]=true;
				}
				else*/
				if ( z[0]+h[0] > z[i]+h[i] ){
					H[i] = z[i] + h[i];
					theta[i] = atan( ((z[0]+h[0]) - (z[i]+h[i])) / w[i] );
					n_eliminated[i]=true;
				}else
					n_eliminated[i]=false;
			}//for

			do {
				loop = false;
				avg = h[0];
				counter = 0;
#pragma loop unroll
				for ( j = 0; j < MOORE_NEIGHBORS; j++)
					if (n_eliminated[j]) {
						avg += H[j];
						counter++;
					}
				if (counter != 0)
					avg = avg / double(counter);
#pragma loop unroll
				for ( j = 0; j < MOORE_NEIGHBORS; j++)
					if (n_eliminated[j] && avg <= H[j]) {
						n_eliminated[j] = false;
						loop = true;
					}
			} while (loop);

			//now collect and save flows
#pragma loop unroll
			for (int i=1;i<MOORE_NEIGHBORS;i++)
				if ( n_eliminated[i] && h[0] > hc*cos(theta[i]) )
				{
					//dd=Pr[i]*(avg - H[i]);
					//d_updateCellValue(current,i+4,dd,x,y);
					//FLOWN-1 return the substate just before the flows.
					//i starts from 1 -> hence it only operates on the flows substates
					d_sbts_updated[d_getIdx(row,col,FLOWN+i-1)]=Pr[i]*(avg - H[i]);
				}
				else
					d_sbts_updated[d_getIdx(row,col,FLOWN+i-1)]=0;



		}//d_sbts_current(THICKNESS)>0

	}//isWithinCABounds
}


__inline__
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
