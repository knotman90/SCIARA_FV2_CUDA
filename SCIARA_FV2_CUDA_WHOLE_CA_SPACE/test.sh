#/bin/bash

echo "THICKNESS"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA/data/2006/output/THICKNESS.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Thickness.stt
echo "TEMPERATURE"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA/data/2006/output/TEMPERATURE.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Temperature.stt
echo "QUOTE"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA/data/2006/output/QUOTE.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_Morphology.stt
echo "SOLIDIFIED"
perl compareMatrix.pl /home/knotman/git/SCIARA_FV2_CUDA/SCIARA_FV2_CUDA/data/2006/output/SOLIFIED_LAVA.out.sst /home/knotman/workspace/C++2/SCIARA_FV2_SERIAL/output/ORIGINAL_000000010000_SolidifiedLavaThickness.stt


