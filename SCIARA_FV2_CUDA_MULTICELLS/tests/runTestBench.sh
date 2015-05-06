#!/bin/bash
#if [ $# -eq 0 ]; then
#    echo "No arguments provided. Needed STEP LIMITS"
#    exit 1
#fi


SERIAL_HOME=/home/knotman/workspace/C++2/SCIARA_FV2_SERIAL
CUDA_HOME=/home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID

SERIAL_EXE=$SERIAL_HOME/Debug/SCIARA_FV2_SERIAL
SERIAL_CONF=$SERIAL_HOME/data201/2006/2006_000000000000.cfg

CUDA_EXE=$CUDA_HOME/Debug/SCIARA_FV2_CUDA_ADAPTIVE_GRID
CUDA_CONF=$CUDA_HOME/data/2006/PARAMETERS.cfg

#LIMIT=$1
#INCREASE=100


RESULT_FILE_SERIAL=SERIAL_OUTPUT.dat
RESULT_FILE_CUDA=CUDA_OUTPUT.dat
touch $RESULT_FILE_SERIAL
echo "TEST RUN FOR $SERIAL_EXE 100 to $LIMIT with step of $INCREASE" > $RESULT_FILE_SERIAL

touch $RESULT_FILE_CUDA
echo "TEST RUN FOR $CUDA_EXE 100 to $LIMIT with step of $INCREASE" > $RESULT_FILE_CUDA

#to generate step use the subsequent haskell instructions (in ghci)
	# let c n = (take (6- (length (n))) $ repeat '0')  ++  n
	# let i= [c ((show n)) | n<-[10,20..100] ]
	# :m + Data.List
	# intercalate " " i
#for ((i=100;i<$LIMIT;i+=$INCREASE)); do
for i in 000500 001000 001500 002000 002500 003000 003500 004000 004500 005000 005500 006000 006500 007000 007500 008000 008500 009000 009500 010000
do
#run the serial version
	#modify the conf file
	sed -i '1s/.*/maximum_steps_(0_for_loop) '$i'/' $SERIAL_CONF
	#run the serial simulation and save results	
	$SERIAL_EXE -c $SERIAL_CONF #>> $RESULT_FILE_SERIAL

	sed -i '1s/.*/nsteps '$i'/' $CUDA_CONF
	$CUDA_EXE -c $CUDA_HOME/data/2006/  #>> $RESULT_FILE_CUDA

	./compareMatrices_steps.sh $i Y Y Y Y > output/comparison/result_comparision_"$i".txt
	
done
