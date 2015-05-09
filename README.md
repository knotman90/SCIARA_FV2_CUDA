# SCIARA FV2 - CUDA IMPLEMENTATION
Sciara FV2 Cellular automaton model for lava flow

Each first level folder of this repository contains a fully working and independent CUDA-nsight (eclipse) project.
IN order to get a working development copy of the model:

1) clone the repo
2) import the project into nsight
3) clean the project (simply delete all the object files from Debug and Release folder if present)
4) set up the cuda setting to match your available hardware (compute capability 3.5 is set by default)
5) Compile.

Execution requires various configuration files as for instance the morphology of the volcano and its lava emission rate during the eruptive event. Configurations are stored in the data folder and can be passed to the executable via command line argument using the -c option and the path of the configuration file (extension .cfg).

For instance ./sciara -c ../data/2006/PARAMETERS.cfg start a simulation based on the eruption of the 2006 at mount Etna in Sicily. PARAMETER.cfg  contains simulation parameters as well as the number of steps for which the simulation have to be carried out.


