################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../MGCNetwork.cc 

OBJS += \
./MGCNetwork.o 

CC_DEPS += \
./MGCNetwork.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../../libraries/CLBllib -I../../../libraries/tnt -I../../../libraries/CNlib2 -I../../../libraries/numlib -I../../../libraries/ISAAC_C++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


