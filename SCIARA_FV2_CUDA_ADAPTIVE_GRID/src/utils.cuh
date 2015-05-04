/*
 * utils.h
 *
 *  Created on: 15/apr/2015
 *      Author: knotman
 */

#ifndef UTILS_H_
#define UTILS_H_
#include<string.h>
#include<iostream>
#include "vent.h"
#include <limits>
#include <vector>

//#####MATH DEFINES####
#define RAD2 (1.41421356237)

//the last enum is just used to have the total number of substates! DO NOT USE IT IN PRODUCTION!
#define VON_NEUMANN_NEIGHBORS	(5)
#define MOORE_NEIGHBORS		(9)

enum SubstatesNames {ALTITUDE=0,THICKNESS,TEMPERATURE,SOLIDIFIED,FLOWN,FLOWO,FLOWE,FLOWS, FLOWNO, FLOWSO, FLOWSE,FLOWNE,NUMBEROFSUBSTATES};
enum AdaptiveGridBoundaries {ROW_START=0,ROW_END,COL_START,COL_END,
	NEW_ROW_START,NEW_ROW_END,NEW_COL_START,NEW_COL_END,ADAPTIVEGRID_SIZE};

/*

		 5 | 1 | 8
		---|---|---
		 2 | 0 | 3
		---|---|---
		 6 | 4 | 7
*/
enum NeighborsIndices {CENTER=0,NORTH,WEST,EAST,SOUTH,NORTH_WEST,SOUTH_WEST,SOUTH_EAST,NORTH_EAST};


inline __host__ __device__
int get_X_LinearIdxVentCoord(int ventIdx){
	return ventIdx*2;
}
inline __host__ __device__
int get_Y_LinearIdxVentCoord(int ventIdx){
	return ventIdx*2+1;
}

/**
 * 2D Coordinates
 */
template <class T>
struct pair{
	T x;
	T y;
};
typedef pair<int> coord;

/** utility function to compute the grid size */
inline
__host__ __device__ int divup(int x, int y) { return x / y + (x % y ? 1 : 0); }


void fatalErrorExit(const char* errmsg){
	fprintf(stderr,"FATAL ERROR: %s\nEXITING",errmsg);
	exit(-1);

}


__host__ __device__ void computeKernelLaunchParameter
(unsigned int threadsBlockX,unsigned int threadsBlockY,
		unsigned int rows,unsigned int cols, dim3 &dimGrid
){
	dimGrid.x = divup(cols,threadsBlockY);
	dimGrid.y = divup(rows,threadsBlockX);



}

/**
 * Compute the number of block on x and y direction necessary to
 * allocate the number of thread described by the adaptive grid bounding box
 * For each axis plus is added to the number of threads.
 * @param dimBlock
 * @param h_d_adaptive_grid (allocated with at least 4 elements describing the adaptive grid (usually unified memory)
 * @param plus
 * @param dimGrid
 */
__host__ __device__ void computeKernelLaunchParameter_plus(dim3 dimBlock,uint* h_d_adaptive_grid,uint plus, dim3 &dimGrid){
	cudaDeviceSynchronize();
	int COLS=h_d_adaptive_grid[COL_END]-h_d_adaptive_grid[COL_START]+plus;
	int ROWS=h_d_adaptive_grid[ROW_END]-h_d_adaptive_grid[ROW_START]+plus;
	dimGrid.x = divup(COLS,dimBlock.x);
	dimGrid.y = divup(ROWS,dimBlock.x);
	//printf("Launch paramrters: BLOCK %i,%i GRID %i,%i\n",dimBlock.x,dimBlock.y,dimGrid.x,dimGrid.y);
	cudaDeviceSynchronize();

}







struct CommandLine{
	//parameters
	bool	_load_config;
	char*	_load_path;
	bool	_verbose;
	bool	_save_config;
	bool	_save_thick;				//new
	char*	_save_path;
	double	_intermediate_fitness;
	bool	_ga_mode;
	bool	_firb_mode;
	int		_np;

	//functions
	void printCmdHelp();
	bool parseArgv(int argc, char* argv[]);
};

void CommandLine::printCmdHelp(){
	printf("USAGE: -c PATH TO CONF FOLDER\n");
}





bool CommandLine::parseArgv(int argc, char* argv[]){
	{
		if(argc>1){

			char config_option_short_str[]	= "-c";
			char config_option_str[]		= "-config";
			char verbose_short_str[]		= "-v";
			char verbose_str[]				= "-verbose";

			char save_thick_short_str[]		= "-st";
			char save_thick_str[]			= "-save_thick";

			char save_step_short_str[]		= "-ss";
			char save_step_str[]			= "-save_step";
			char firb_mode_str[]			= "-firb_mode";
			char np_str[]					= "-np";

			int i = 2;
			while(argv[i])
			{
				if (argv[i][0] == '-')
					if (
							strcmp(argv[i], config_option_short_str)	&&
							strcmp(argv[i], config_option_str)			&&
							strcmp(argv[i], verbose_short_str)			&&
							strcmp(argv[i], verbose_str)				&&
							strcmp(argv[i], save_thick_short_str)		&&
							strcmp(argv[i], save_thick_str)			&&
							strcmp(argv[i], save_step_short_str)		&&
							strcmp(argv[i], save_step_str)				&&
							strcmp(argv[i], firb_mode_str)				&&
							strcmp(argv[i], np_str)
					){

						printCmdHelp();
						fatalErrorExit("BAD ARGUMENTS");

					}//if

				i = i + 1;
			}

			i = 1;
			while(argv[i])
			{
				if (!strcmp(argv[i], config_option_short_str) || !strcmp(argv[i], config_option_str))
				{
					_load_config = true;
					_load_path = argv[i+1];
				}
				if (!strcmp(argv[i], verbose_short_str) || !strcmp(argv[i], verbose_str))
					_verbose = true;

				if (!strcmp(argv[i], firb_mode_str))
					_firb_mode = true;

				if (!strcmp(argv[i], np_str))
					_np = atoi(argv[i+1]);

				i = i + 1;
			}
		}
	}
	return true;
}


/**
 * Adaptivegrid array should be allocated (ADAPTIVEGRID_SIZE size)
 * Construct the minimum bounding box for that containts the vents
 * @param adaptiveGrid
 * @param vent
 */
void initializeadaptiveGrid(uint* adaptiveGrid,vector<TVent> vent){
	adaptiveGrid[NEW_ROW_START]	=adaptiveGrid[ROW_START] = UINT_MAX;
	adaptiveGrid[NEW_COL_START]	=adaptiveGrid[COL_START] = UINT_MAX;

	adaptiveGrid[NEW_ROW_END]		=adaptiveGrid[ROW_END]   = 0;
	adaptiveGrid[NEW_COL_END]		=adaptiveGrid[COL_END]   = 0;

	for(auto v : vent){
		adaptiveGrid[NEW_ROW_START]	=adaptiveGrid[ROW_START] = min(adaptiveGrid[ROW_START],v.y());
		adaptiveGrid[NEW_ROW_END]	=adaptiveGrid[ROW_END]   = max(adaptiveGrid[ROW_END],v.y());

		adaptiveGrid[NEW_COL_START]	=adaptiveGrid[COL_START] = min(adaptiveGrid[COL_START],v.x());
		adaptiveGrid[NEW_COL_END]	=adaptiveGrid[COL_END]   = max(adaptiveGrid[COL_END],v.x());

	}
}

#endif /* UTILS_H_ */
