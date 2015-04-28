#/bin/bash

echo "--SCIARA MATRIX TESTS--"
echo "Takes max 4 parameters: if the k'th command line argument is present then the k'th test is verbose"
echo "THICKNESS"
perl compareMatrix.pl  /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/THICKNESS.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Thickness.stt $1
echo "TEMPERATURE"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/TEMPERATURE.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Temperature.stt $2
echo "QUOTE"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/QUOTE.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Morphology.stt $3
echo "SOLIDIFIED"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA_ADAPTIVE_GRID/data/2006/output/SOLIFIED_LAVA.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_SolidifiedLavaThickness.stt $4


