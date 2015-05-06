################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/matrix.cpp \
../src/vent.cpp 

C_SRCS += \
../src/scratch.c 

CU_SRCS += \
../src/cuda_error_check.cu \
../src/main.cu 

CU_DEPS += \
./src/cuda_error_check.d \
./src/main.d 

OBJS += \
./src/cuda_error_check.o \
./src/main.o \
./src/matrix.o \
./src/scratch.o \
./src/vent.o 

C_DEPS += \
./src/scratch.d 

CPP_DEPS += \
./src/matrix.d \
./src/vent.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 -gencode arch=compute_35,code=sm_35  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 -gencode arch=compute_35,code=sm_35  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 -gencode arch=compute_35,code=sm_35  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.0/bin/nvcc -G -g -O0 -Xcompiler -fPIC -std=c++11 --compile  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


