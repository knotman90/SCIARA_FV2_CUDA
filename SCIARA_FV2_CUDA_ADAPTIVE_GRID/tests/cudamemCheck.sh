#Author: Davide Spataro

#/bin/bash

if [ $# -eq 0 ]; then
    echo "No arguments provided. Please specify r o d for release or debug respectively"
    exit 1
fi

MEMCHECK_BIN=/usr/local/cuda-7.0/bin/cuda-memcheck
MEMCHECK_DEFAULT_OPTIONS="--language c"
MEMCHECK_OPTIONS_DEFAULT="--tool memcheck "
MEMCHECK_OPTIONS_LEAK="--leak-check full"
MEMCHECK_OPTIONS_SHAREDRACECHECK="--tool racecheck"
MEMCHECK_OPTIONS_INITCHECK="--tool initcheck"
MEMCHECK_OPTIONS_SYNCHCHECK="--tool synccheck"

APP_HOME=/home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID
APP_BIN=""
APP_OPTIONS="-c $APP_HOME/data/2006/"

if [ $1 == "d" ]; then
	APP_BIN=Debug
elif [ $1 == "r" ]; then
	APP_BIN=Release
else
	echo "Please specify r or s as first argument"
    exit 1
fi

APP_BIN=$APP_HOME/$APP_BIN/SCIARA_FV2_CUDA_ADAPTIVE_GRID
#-----------------------DEFAULT-------------------------------------------
	#The tool can precisely detect and report out of bounds and misaligned memory accesses to
	#global, local, shared and global atomic instructions in CUDA applications
echo ""
echo "Performing DEFAULT MEMCHECK on application binary"
$MEMCHECK_BIN $MEMCHECK_OPTIONS_DEFAULT $MEMCHECK_DEFAULT_OPTIONS $APP_BIN $APP_OPTIONS
echo ""
#-----------------------LEAK-------------------------------------------
	#The memcheck tool can detect leaks of allocated memory.
	#Memory leaks are device side allocations that have not 
	#been freed by the time the context is destroyed.
echo ""
echo "Performing  LEAK on application binary "
$MEMCHECK_BIN $MEMCHECK_OPTIONS_LEAK $MEMCHECK_DEFAULT_OPTIONS $APP_BIN $APP_OPTIONS
echo ""
echo ""
#-----------------------RACE CHECK-------------------------------------------
	#The racecheck tool is a run time shared memory data access hazard detector. 
	#The primary use of this tool is to help identify memory access race 
	#conditions in CUDA applications that use shared memory.
echo ""
echo "Performing  RACECHECK on application binary "
$MEMCHECK_BIN $MEMCHECK_OPTIONS_SHAREDRACECHECK $MEMCHECK_DEFAULT_OPTIONS $APP_BIN $APP_OPTIONS
echo ""
echo ""
#-----------------------INIT CHECK-------------------------------------------
	#The initcheck tool is a run time uninitialized device global memory access detector.
echo ""
echo "Performing  INITCHECK on application binary "
$MEMCHECK_BIN $MEMCHECK_OPTIONS_INITCHECK $MEMCHECK_DEFAULT_OPTIONS $APP_BIN $APP_OPTIONS
echo ""
echo ""

#-----------------------SYNCH CHECK-------------------------------------------
	#The synccheck tool is a runtime tool that can identify whether a CUDA
	#application is attempting to call __syncthreads() from divergent code.
echo ""
echo "Performing  SYNCH on application binary "
$MEMCHECK_BIN $MEMCHECK_OPTIONS_SYNCHCHECK $MEMCHECK_DEFAULT_OPTIONS $APP_BIN $APP_OPTIONS
echo ""
echo ""



