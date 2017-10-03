################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../core/filter.cpp \
../core/flipImage.cpp \
../core/hough.cpp \
../core/houghUtilities.cpp 

OBJS += \
./core/filter.o \
./core/flipImage.o \
./core/hough.o \
./core/houghUtilities.o 

CPP_DEPS += \
./core/filter.d \
./core/flipImage.d \
./core/hough.d \
./core/houghUtilities.d 


# Each subdirectory must supply rules for building sources it contributes
core/%.o: ../core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


