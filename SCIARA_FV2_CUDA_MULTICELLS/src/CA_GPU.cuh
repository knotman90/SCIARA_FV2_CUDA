/*
 * CA_GPU.cuh
 *
 *  Created on: 15/apr/2015
 *      Author: knotman
 */

#ifndef CA_GPU_CUH_
#define CA_GPU_CUH_
#include <math.h>
#define DEB (0)

#define BDIM_X (16)
#define BDIM_Y (16)

class CA_GPU{
	friend class CA_HOST;
public:
	//ADAPTIVE GRID
	uint* h_d_adaptive_grid; //managed cuda Unified address
private:


	//#############################  CA VARIABLES (DEVICE)  ##########################
	double d_sim_elapsed_time; //tempo trascorso dall'inizio della simulazione [s]

	//#############################  CA PARAMETERS (DEVICE)  ########################
	//																				#
	//																				#
	int d_nSteps;//																	#
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
	double* d_sbts_updated; //linearized substates
	double* d_sbts_current; //linearized substates

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
	__device__ int getMooreNeighIdx_ROW(int row, int neighIdx);

	__inline__
	__device__ int getMooreNeighIdx_COL(int col, int neighIdx);

	__inline__
	__device__ int d_getIdx(int x, int y, int substate)	{
		return ( (d_NUMCELLS * substate)  +   (x * d_NC) + y );
	}

	__inline
	__device__ bool isWithinCABounds(int row, int col);

	__inline
	__device__ bool isWithinCABounds_WholeSpace(int row, int col);

	__inline
	__device__ bool isWithinCABounds_AG(int row, int col);

	__inline
	__device__ bool isWithinCABounds_AG_FAT(int row, int col);


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

	__inline__
	__device__ double PowerLaw(double k1, double k2, double T){
		return exp10(k1 + k2*T);
	}

	__device__ void distribuiteFlows();

	__inline__
	__device__ double d_computeNewTemperature(double sommah, double sommath);


	__device__ void swapMatrices();


	__inline__
	__device__ void adjustAdaptiveGrid(int row, int col);

	__inline__
	__device__ void copyNewToCurrentAdaptiveGrid();


};//end definition of CA_GPU


__device__ void CA_GPU::swapMatrices(){
#pragma unroll
	for(int cycle=0;cycle<CPT_COL;cycle++){

		int row = blockIdx.y * blockDim.y + threadIdx.y+h_d_adaptive_grid[ROW_START];
		int col = blockIdx.x * CPT_COL*(blockDim.x) + (cycle*blockDim.x) + threadIdx.x+h_d_adaptive_grid[COL_START];
		int idx;
		if(isWithinCABounds(row,col)){

			//Thickness
			idx=d_getIdx(row,col,ALTITUDE);
			d_sbts_updated[idx] = d_sbts_current[idx];
			//Thickness
			idx=d_getIdx(row,col,THICKNESS);
			d_sbts_updated[idx] = d_sbts_current[idx];
			//Temperature
			idx=d_getIdx(row,col,TEMPERATURE);
			d_sbts_updated[idx] = d_sbts_current[idx];
			//Solidified
			idx=d_getIdx(row,col,SOLIDIFIED);
			d_sbts_updated[idx] = d_sbts_current[idx];

			//clean up flows
			for(int i=FLOWN;i<=FLOWNE;i++){
				idx=d_getIdx(row,col,i);
				d_sbts_current[idx]=0.0;
				d_sbts_updated[idx]=0.0;
			}
		}
	}

	return;
}
#define EPISILON (0.0)
//TODO unroll sums caching the flows values
__device__ void CA_GPU::distribuiteFlows(){
#pragma unroll
	for(int cycle=0;cycle<CPT_COL;cycle++){

		//get cell coordinates
		int row = blockIdx.y * blockDim.y + threadIdx.y+h_d_adaptive_grid[ROW_START]-1;
		int col = blockIdx.x * CPT_COL*(blockDim.x) + (cycle*blockDim.x) + threadIdx.x+h_d_adaptive_grid[COL_START]-1;

		if(isWithinCABounds_AG_FAT(row,col)){
			//newtick = subtract(sum of outflows) + sum(inflows)

			double sommah=d_sbts_updated[d_getIdx(row,col,THICKNESS)];
			double new_temp = d_sbts_updated[d_getIdx(row,col,TEMPERATURE)];
			double sommath;

			//flown(x,y) has to be added comes from the sud cell
			//flows(neighbor1) has to be substracted because is the amount of lava that the current cell
			//transfer to the upper cell
			//same hold for the other pair of indices

			//new_thick= -flowToNeigh + flowReceived
			if(isWithinCABounds(row,col)){
				sommah-=
						d_sbts_current[d_getIdx(row,col,FLOWN)]+
						d_sbts_current[d_getIdx(row,col,FLOWS)]+
						d_sbts_current[d_getIdx(row,col,FLOWE)]+
						d_sbts_current[d_getIdx(row,col,FLOWO)]+
						d_sbts_current[d_getIdx(row,col,FLOWNO)]+
						d_sbts_current[d_getIdx(row,col,FLOWNE)]+
						d_sbts_current[d_getIdx(row,col,FLOWSO)]+
						d_sbts_current[d_getIdx(row,col,FLOWSE)];
			}


			sommath= sommah * new_temp;
			sommah+=
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH),getMooreNeighIdx_COL(col,NORTH),FLOWS)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH),getMooreNeighIdx_COL(col,SOUTH),FLOWN)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,WEST),getMooreNeighIdx_COL(col,WEST),FLOWE)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,EAST),getMooreNeighIdx_COL(col,EAST),FLOWO)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_WEST),getMooreNeighIdx_COL(col,SOUTH_WEST),FLOWNE)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_WEST),getMooreNeighIdx_COL(col,NORTH_WEST),FLOWSE)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_EAST),getMooreNeighIdx_COL(col,SOUTH_EAST),FLOWNO)]+
					d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_EAST),getMooreNeighIdx_COL(col,NORTH_EAST),FLOWSO)];

			sommath+=
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH),getMooreNeighIdx_COL(col,NORTH),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH),getMooreNeighIdx_COL(col,NORTH),FLOWS)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH),getMooreNeighIdx_COL(col,SOUTH),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH),getMooreNeighIdx_COL(col,SOUTH),FLOWN)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,WEST),getMooreNeighIdx_COL(col,WEST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,WEST),getMooreNeighIdx_COL(col,WEST),FLOWE)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,EAST),getMooreNeighIdx_COL(col,EAST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,EAST),getMooreNeighIdx_COL(col,EAST),FLOWO)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_WEST),getMooreNeighIdx_COL(col,SOUTH_WEST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_WEST),getMooreNeighIdx_COL(col,SOUTH_WEST),FLOWNE)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_WEST),getMooreNeighIdx_COL(col,NORTH_WEST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_WEST),getMooreNeighIdx_COL(col,NORTH_WEST),FLOWSE)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_EAST),getMooreNeighIdx_COL(col,SOUTH_EAST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_EAST),getMooreNeighIdx_COL(col,SOUTH_EAST),FLOWNO)]
					)+
					(d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_EAST),getMooreNeighIdx_COL(col,NORTH_EAST),TEMPERATURE)]*
							d_sbts_current[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_EAST),getMooreNeighIdx_COL(col,NORTH_EAST),FLOWSO)]);


			if(sommah > EPISILON){
				//update thickness and temperature
				double new_temp = d_computeNewTemperature(sommah,sommath);
				d_sbts_current[d_getIdx(row,col,TEMPERATURE)]= new_temp;
				d_sbts_current[d_getIdx(row,col,THICKNESS)]	 = sommah;

				//quote increment due to solidification
				if(new_temp < d_PTsol){
					if(DEB)printf("Solidified %i,%i %.5f %.5f %.5f\n",row,col,d_sbts_current[d_getIdx(row,col,THICKNESS)],new_temp,d_PTsol);
					double newQuote = d_sbts_updated[d_getIdx(row,col,ALTITUDE)]+d_sbts_current[d_getIdx(row,col,THICKNESS)];
					double newSolid = d_sbts_updated[d_getIdx(row,col,SOLIDIFIED)]+d_sbts_current[d_getIdx(row,col,THICKNESS)];
					d_sbts_current[d_getIdx(row,col,SOLIDIFIED)] = newSolid;
					d_sbts_current[d_getIdx(row,col,ALTITUDE)] = newQuote;
					d_sbts_current[d_getIdx(row,col,THICKNESS)] = 0;
					d_sbts_current[d_getIdx(row,col,TEMPERATURE)] = d_PTsol;
				}else{
					//there is lava and is not solidified -> activate this cell!
					adjustAdaptiveGrid(row,col);
				}
			}
		}
	}//for
}


__inline__
__device__ double CA_GPU::d_computeNewTemperature(double sommah, double sommath){
	double new_temp = sommath/sommah;
	double aus = 1.0 + (3 *  pow(new_temp, 3.0) * d_Pepsilon * d_Psigma * d_Pclock * d_Pcool) / (d_Prho * d_Pcv * sommah * d_Pac);
	new_temp = new_temp /pow(aus,1.0/3.0);
	return new_temp;
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
__device__ int CA_GPU::getMooreNeighIdx_ROW(int row, int neighIdx){
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
__device__ int CA_GPU::getMooreNeighIdx_COL(int col, int neighIdx){
	if(neighIdx==0 || neighIdx==1 || neighIdx == 4)
		return col;

	if(neighIdx==2 || neighIdx==5 || neighIdx == 6)
		return col-1;

	//implicit else if(neighIdx==3 || neighIdx==7 || neighIdx == 8)
	return col+1;

}


__inline
__device__ bool CA_GPU::isWithinCABounds(int row, int col){

	//return isWithinCABounds_WholeSpace(row,col);
	return isWithinCABounds_AG(row,col);
}



__inline
__device__ bool CA_GPU::isWithinCABounds_WholeSpace(int row, int col){

	return
			(row>=0 && row <= (d_NR-1)) &&
			col>=0 && col <= (d_NC-1);
}

/**
 * Chack if a threads is in the boundaries
 * described by the adaptive grid boundaries
 * ADDING A BORDER OF 1 cell
 * h_d_adaptive_grid
 * @param row
 * @param col
 * @param
 * @return
 */
__inline
__device__ bool CA_GPU::isWithinCABounds_AG_FAT(int row, int col ){

	return
			(row>=h_d_adaptive_grid[ROW_START]-1 && row <= h_d_adaptive_grid[ROW_END]+1) &&
			col>=h_d_adaptive_grid[COL_START]-1 && col <= h_d_adaptive_grid[COL_END]+1;
}

/**
 * Chack if a threads is in the boundaries
 * described by the adaptive grid boundaries
 * h_d_adaptive_grid
 * @param row
 * @param col
 * @param
 * @return
 */
__inline
__device__ bool CA_GPU::isWithinCABounds_AG(int row, int col ){

	return
			(row>=h_d_adaptive_grid[ROW_START] && row <= h_d_adaptive_grid[ROW_END]) &&
			col>=h_d_adaptive_grid[COL_START] && col <= h_d_adaptive_grid[COL_END];
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
		int x = coordVents[get_X_LinearIdxVentCoord(vent)];
		int y = coordVents[get_Y_LinearIdxVentCoord(vent)];

		int cellIdx=d_getIdx(x,y,THICKNESS);
		d_sbts_updated[cellIdx] += emitted_lava;
		d_sbts_current[cellIdx] = d_sbts_updated[cellIdx];

		cellIdx=d_getIdx(x,y,TEMPERATURE);
		d_sbts_updated[cellIdx] = d_PTvent;
		d_sbts_current[cellIdx] = d_PTvent;

		//if(DEB) printf("emitting lava (%d,%d), %f\n",x,y,emitted_lava);
	}

	if(threadIdx.x==0 && threadIdx.y==0){
		d_sim_elapsed_time+=d_Pclock;
	}

}


__inline__
__device__ void CA_GPU::cellTemperatureInitialize(){
#pragma unroll
	for(int cycle=0;cycle<CPT_COL;cycle++){

		//get cell coordinates
		int row = blockIdx.y * blockDim.y + threadIdx.y+h_d_adaptive_grid[ROW_START];
		int col = blockIdx.x * CPT_COL*(blockDim.x) + (cycle*blockDim.x) + threadIdx.x+h_d_adaptive_grid[COL_START];

		if(isWithinCABounds(row,col)){
			//current cell temperature
			double thickness = d_sbts_updated[d_getIdx(row,col,THICKNESS)];

			if(thickness > 0){
				double temp = d_sbts_updated[d_getIdx(row,col,TEMPERATURE)];
				d_sbts_current[d_getIdx(row,col,TEMPERATURE)] = thickness*temp;
			}
			return;
		}

	}
}

__device__ void CA_GPU::empiricalFlows(){
#pragma unroll
	for(int cycle=0;cycle<CPT_COL;cycle++){
	__syncthreads();

		int row = blockIdx.y * blockDim.y + threadIdx.y+h_d_adaptive_grid[ROW_START];
		int col = blockIdx.x * CPT_COL*(blockDim.x) + (cycle*blockDim.x) + threadIdx.x+h_d_adaptive_grid[COL_START];

		/*
		 * Shared memory has size dimBlock.x+2 * dimBlock.y+2.
		 * The indices of the locations that belongs to threads of this block are from the first upper
		 * left c (1,1) to (dimBlock.x,dimBlock.y). This means that each cells can safely and easily fill
		 * its indices (just translated by one row and column) and border have to be managed separately.
									    SHARED MEMORY LOGIC
											+--------+
											|uffffffy|
											|kcccccch|
											|kcccccch|
											|kcccccch|
											|kcccccch|
											|kcccccch|
											|kcccccch|
											|tllllllr|
											+--------+
					X are ghost cells (values of threads of neighbors blocks)
					c are cells inside the block
		 */
		/**The row index of the current thread in the shared memory matrix*/
		int row_s=threadIdx.y+1;
		/**The column index of the current thread in the shared memory matrix*/
		int col_s=threadIdx.x+1;
		//### SHARED INITIALIZATION ###
		__shared__ double s_altitude [BDIM_Y+2][BDIM_X/CPT_COL+2];
		__shared__ double s_thickness[BDIM_Y+2][BDIM_X/CPT_COL+2];

		//-------------------FIRST ghost row (f in the picture)
		if(threadIdx.y==0){
			s_altitude	[0][col_s]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH),getMooreNeighIdx_COL(col,NORTH),ALTITUDE)] ;
			s_thickness	[0][col_s]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH),getMooreNeighIdx_COL(col,NORTH),THICKNESS)] ;

			//upperleft corner also initialized by the thread (0,0) (u in picture)
			if(threadIdx.x==0){
				s_altitude	[0][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_WEST),getMooreNeighIdx_COL(col,NORTH_WEST),ALTITUDE)];
				s_thickness	[0][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_WEST),getMooreNeighIdx_COL(col,NORTH_WEST),THICKNESS)];
			}

			//uppertigh corner
			if(threadIdx.x==blockDim.x-1){//upper right corner (y in picture)
				s_altitude	[0][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_EAST),getMooreNeighIdx_COL(col,NORTH_EAST),ALTITUDE)];
				s_thickness	[0][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,NORTH_EAST),getMooreNeighIdx_COL(col,NORTH_EAST),THICKNESS)];
			}
		}

		//-------------------LAST ghost row (l in the picture)
		if(threadIdx.y==blockDim.y-1){
			s_altitude	[blockDim.y+1][col_s]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH),getMooreNeighIdx_COL(col,SOUTH),ALTITUDE)] ;
			s_thickness	[blockDim.y+1][col_s]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH),getMooreNeighIdx_COL(col,SOUTH),THICKNESS)] ;

			if(threadIdx.x==0){//bottom left corner (t in picture)
				s_altitude	[blockDim.y+1][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_WEST),getMooreNeighIdx_COL(col,SOUTH_WEST),ALTITUDE)] ;
				s_thickness	[blockDim.y+1][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_WEST),getMooreNeighIdx_COL(col,SOUTH_WEST),THICKNESS)] ;
			}
			if(threadIdx.x==blockDim.x-1){//bottom right corner (r in picture)
				s_altitude	[blockDim.y+1][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_EAST),getMooreNeighIdx_COL(col,SOUTH_EAST),ALTITUDE)] ;
				s_thickness	[blockDim.y+1][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,SOUTH_EAST),getMooreNeighIdx_COL(col,SOUTH_EAST),THICKNESS)] ;
			}
		}

		//-------------------LEFT ghost COLUMN  (k in the picture)
		if(threadIdx.x==0){
			s_altitude	[row_s][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,WEST),getMooreNeighIdx_COL(col,WEST),ALTITUDE)] ;
			s_thickness	[row_s][0]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,WEST),getMooreNeighIdx_COL(col,WEST),THICKNESS)] ;
		}
		//-------------------LEFT ghost COLUMN  (h in the picture)
		if(threadIdx.x==blockDim.x-1){
			s_altitude	[row_s][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,EAST),getMooreNeighIdx_COL(col,EAST),ALTITUDE)] ;
			s_thickness	[row_s][blockDim.x+1]	= d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,EAST),getMooreNeighIdx_COL(col,EAST),THICKNESS)] ;
		}



		//cells inside the block done by all threads of the block! (c in picture)
		s_altitude	[row_s][col_s]	= d_sbts_updated[d_getIdx(row,col,ALTITUDE)];
		s_thickness	[row_s][col_s]	= d_sbts_updated[d_getIdx(row,col,THICKNESS)];

		//### SHARED SYNCH ###
		__syncthreads(); //shared fence -> data has to be correctly initialized by all threads



		if(isWithinCABounds(row,col)){
			if (s_thickness[row_s][col_s] > 0) {

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
				_Pr = PowerLaw(d_a, d_b, d_sbts_updated[d_getIdx(row,col,TEMPERATURE)]);
				hc 	= PowerLaw(d_c, d_d, d_sbts_updated[d_getIdx(row,col,TEMPERATURE)]);


#pragma unroll
				for ( i = 0; i < MOORE_NEIGHBORS; i++) {

					//h[i] = d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,i),getMooreNeighIdx_COL(col,i),THICKNESS)] ;
					h[i] = s_thickness[getMooreNeighIdx_ROW(row_s,i)][getMooreNeighIdx_COL(col_s,i)];
					H[i]  = theta[i] = 0;
					w[i] = _w;
					Pr[i] = _Pr;

					if (i < VON_NEUMANN_NEIGHBORS){
						z[i] = s_altitude[getMooreNeighIdx_ROW(row_s,i)][getMooreNeighIdx_COL(col_s,i)];
						//z[i] = d_sbts_updated[d_getIdx(getMooreNeighIdx_ROW(row,i),getMooreNeighIdx_COL(col,i),ALTITUDE)] ;
					}else{
						//val - (val-valNeigh_I)/rad2
						z[i]= s_altitude[row_s][col_s] -
								(
										s_altitude[row_s][col_s] -
										s_altitude[getMooreNeighIdx_ROW(row_s,i)][getMooreNeighIdx_COL(col_s,i)]

								)/RAD2;

					}//else
				}//for

				H[0] = z[0];
				n_eliminated[0]=true;

#pragma unroll
				for ( i = 1; i < MOORE_NEIGHBORS; i++){

					if ( z[0]+h[0] > z[i]+h[i] ){
						H[i] = z[i] + h[i];
						theta[i] = atan( ((z[0]+h[0]) - (z[i]+h[i])) / w[i] );
						n_eliminated[i]=true;
					}else
						n_eliminated[i]=false;

					//if(DEB) printf("(%i,%i),z%i=%.5f, eliminated=%i, theta=%.5f, H=%.5f, h=%.5f, w=%.5f\n",row,col,i,z[i],n_eliminated[i],theta[i],H[i],h[i],w[i]);
				}//for

				do {
					loop = false;
					avg = h[0];
					counter = 0;
#pragma unroll
					for ( j = 0; j < MOORE_NEIGHBORS; j++)
						if (n_eliminated[j]) {
							avg += H[j];
							counter++;
						}
					if (counter != 0)
						avg = avg / double(counter);
#pragma unroll
					for ( j = 0; j < MOORE_NEIGHBORS; j++)
						if (n_eliminated[j] && avg <= H[j]) {
							n_eliminated[j] = false;
							loop = true;
						}
				} while (loop);

				//now collect and save flows
#pragma unroll
				for (int i=1;i<MOORE_NEIGHBORS;i++)
					if ( n_eliminated[i] && h[0] > hc*cos(theta[i]) )
					{
						//if(DEB) printf("Flusso verso (%d,%d),%.9f,%.9f,%.9f, %.5f flow to %i, temp=%.5f\n",row,col,avg,Pr[i],H[i],Pr[i]*(avg - H[i]),i,d_sbts_updated[d_getIdx(row,col,TEMPERATURE)]);
						//dd=Pr[i]*(avg - H[i]);
						//d_updateCellValue(current,i+4,dd,x,y);
						//FLOWN-1 return the substate just before the flows.
						//i starts from 1 -> hence it only operates on the flows substates
						d_sbts_current[d_getIdx(row,col,FLOWN+i-1)]= Pr[i]*(avg - H[i]);
					}
					else
						d_sbts_current[d_getIdx(row,col,FLOWN+i-1)]= 0.0;

			}//d_sbts_updated(THICKNESS)>0

		}//isWithinCABounds
	}//for cycle
}



__inline__
__device__ void CA_GPU::adjustAdaptiveGrid(int row, int col) {
	atomicMax(&h_d_adaptive_grid[NEW_COL_END],col);
	atomicMin(&h_d_adaptive_grid[NEW_COL_START],col);

	atomicMax(&h_d_adaptive_grid[NEW_ROW_END],row);
	atomicMin(&h_d_adaptive_grid[NEW_ROW_START],row);

}

/**
 * Should always be launched with ADAPTIVEGRID_SIZE/2 threads!
 */
__inline__
__device__ void CA_GPU::copyNewToCurrentAdaptiveGrid() {
	int ADAPTIVEGRID_SIZE_2=ADAPTIVEGRID_SIZE/2;
	if(blockIdx.x==0 && blockIdx.y==0){
		if(threadIdx.x <=ADAPTIVEGRID_SIZE_2 && threadIdx.y==0){
			h_d_adaptive_grid[threadIdx.x] = h_d_adaptive_grid[threadIdx.x+ADAPTIVEGRID_SIZE_2];
		}
	}

}


__inline__
__device__ void CA_GPU::printSubstate(int substate){
	int row = blockIdx.y * blockDim.y + threadIdx.y+h_d_adaptive_grid[ROW_START];
	int col = blockIdx.x * blockDim.x + threadIdx.x+h_d_adaptive_grid[COL_START];

	if(isWithinCABounds(row,col)){
		double val = d_sbts_updated[d_getIdx(row,col,substate)];
		if(val>0)
			printf("(%d,%d,%d,%d,%d,%d),%.9f \n",row,col,d_getIdx(row,col,substate),d_NUMCELLS,d_NR,d_NC,val);
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
