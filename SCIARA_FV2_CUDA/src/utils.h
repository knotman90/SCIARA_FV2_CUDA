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


struct CommandLine{
	//parameters
	int		_argc;
	char**	_argv;
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
	bool parseArgv(int argc, char** argv);
};

void CommandLine::printCmdHelp(){
	printf("USAGE: -c PATH TO CONF FOLDER\n");
}

bool CommandLine::parseArgv(int argc, char** argv){
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
			char ga_mode_str[]				= "-ga_mode";
			char firb_mode_str[]			= "-firb_mode";
			char np_str[]					= "-np";

			int i = 1;
			while(_argv[i])
			{
				if (_argv[i][0] == '-')
					if (
							strcmp(_argv[i], config_option_short_str)	&&
							strcmp(_argv[i], config_option_str)			&&
							strcmp(_argv[i], verbose_short_str)			&&
							strcmp(_argv[i], verbose_str)				&&
							strcmp(_argv[i], save_thick_short_str)		&&
							strcmp(_argv[i], save_thick_str)			&&
							strcmp(_argv[i], save_step_short_str)		&&
							strcmp(_argv[i], save_step_str)				&&
							strcmp(_argv[i], firb_mode_str)				&&
							strcmp(_argv[i], np_str)
					){

						printf("BAD ARGUMENTS\n");
						printCmdHelp();
						printf("EXITING\n");
						exit(-1);
					}//if

				i = i + 1;
			}

			i = 1;
			while(_argv[i])
			{
				if (!strcmp(_argv[i], config_option_short_str) || !strcmp(_argv[i], config_option_str))
				{
					_load_config = true;
					_load_path = _argv[i+1];
				}
				if (!strcmp(_argv[i], verbose_short_str) || !strcmp(_argv[i], verbose_str))
					_verbose = true;

				if (!strcmp(_argv[i], firb_mode_str))
					_firb_mode = true;

				if (!strcmp(_argv[i], np_str))
					_np = atoi(_argv[i+1]);

				i = i + 1;
			}
		}
	}
}

#endif /* UTILS_H_ */
