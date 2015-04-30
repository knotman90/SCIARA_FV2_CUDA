#/bin/bash

echo ""
echo "---------------------------------------------"
echo "         --SCIARA MATRIX TESTS--"
echo "USAGE: steps (6 digits) and 4 arguments for"
echo "the verbosity of the corrensponding tests)"
echo "./tests 000010 Y Y will perform tests on" 
echo "matrices serial of steps 10 and verbose" 
echo "for thickness and temperature"
echo "--------------------------------------------"
echo ""

STEPS=$1;
SIZE_STEPS="${#STEPS}"
if [[ (-z "$STEPS") ]]
	 then 
	 	echo "NUMBER OF STEPS IS REQUIRED (6 digits)";	exit
fi
if [[ $SIZE_STEPS -ne 6 ]]
	 then 
	 	echo "NUMBER OF STEPS LENGTH HAS TO BE 6 digits";	exit
fi





echo "		THICKNESS"
perl ../compareMatrix.pl  /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/THICKNESS.out.sst \
						"./output/ORIGINAL_000000"$STEPS"_Thickness.stt" $2

echo "		TEMPERATURE"
perl ../compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/TEMPERATURE.out.sst \
						"./output/ORIGINAL_000000"$STEPS"_Temperature.stt" $3

echo "		QUOTE"
perl ../compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/QUOTE.out.sst \
						"./output/ORIGINAL_000000"$STEPS"_Morphology.stt" $4
echo "		SOLIDIFIED"
perl ../compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/SOLIFIED_LAVA.out.sst \
						"./output/ORIGINAL_000000"$STEPS"_SolidifiedLavaThickness.stt" $5
