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

//#####MATH DEFINES####
#define RAD2 (1.41421356237)

//the last enum is just used to have the total number of substates! DO NOT USE IT IN PRODUCTION!
#define VON_NEUMANN_NEIGHBORS	(5)
#define MOORE_NEIGHBORS		(9)

enum SubstatesNames {ALTITUDE=0,THICKNESS,TEMPERATURE,SOLIDIFIED,FLOWN,FLOWO,FLOWE,FLOWS, FLOWNO, FLOWSO, FLOWSE,FLOWNE,NUMBEROFSUBSTATES};


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
__host__ int divup(int x, int y) { return x / y + (x % y ? 1 : 0); }


void fatalErrorExit(const char* errmsg){
	fprintf(stderr,"FATAL ERROR: %s\nEXITING",errmsg);
	exit(-1);

}


void computeKernelLaunchParameter(unsigned int threadsBlockX,unsigned int threadsBlockY,
								  unsigned int cellsX,unsigned int cellsY, dim3 &dimGrid
								){

dimGrid.x = divup(cellsX,threadsBlockX);
dimGrid.y = divup(cellsY,threadsBlockY);


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

#endif /* UTILS_H_ */
