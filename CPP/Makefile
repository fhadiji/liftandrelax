LIB_DAI := /home/stream/lib/libDAI/libDAI-git
PYTHON := /usr/include/python2.6
BOOST_INCLUDE := /usr/include/boost
BOOST_LIB := /usr/lib

RM := rm -rf

OBJS := \
./src/main/AnyPositionCnfCompress.o \
./src/main/AnyPositionCompress.o \
./src/main/CBP.o \
./src/main/CFactorGraph.o \
./src/main/Compress.o \
./src/main/PositionCompress.o \

CPP_DEPS := \
./src/main/AnyPositionCnfCompress.d \
./src/main/AnyPositionCompress.d \
./src/main/CBP.d \
./src/main/CFactorGraph.d \
./src/main/Compress.d \
./src/main/PositionCompress.d \

# All Target
all: libSTREAM.so libSTREAMWrapper.so

# rule for Python wrapper object
src/main/STREAMWrapper.o : 
	@echo 'Building file: src/main/libSTREAMWrapper.cpp'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DDAI_WITH_BP -DDEMO_DOCSTRING_SHOW_ALL=true -I"./include" -I$(LIB_DAI)"/include" -I$(BOOST_INCLUDE) -I$(PYTHON) -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "src/main/STREAMWrapper.cpp"
	@echo 'Finished building: $<'
	@echo ' '

# rule for library objects
src/main/%.o: src/main/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"./include" -I$(LIB_DAI)"/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

libSTREAM.so: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -shared -o"libSTREAM.so" $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

libSTREAMWrapper.so : src/main/STREAMWrapper.o
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -shared -o"libSTREAMWrapper.so" src/main/STREAMWrapper.o libSTREAM.so $(BOOST_LIB)/libboost_python.so -lutil -lpthread -ldl -L$(LIB_DAI)/lib -ldai 
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) $(CPP_DEPS) libSTREAM.so ./src/main/STREAMWrapper.d ./src/main/STREAMWrapper.o libSTREAMWrapper.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


