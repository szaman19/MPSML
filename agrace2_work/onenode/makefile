FLAGS= -std=c++11 -g
PETSC_ARCH = arch_linux_gnu_c_debug
MYLIBS = ${SLEPC_LIB}/libslepc.so ${PETSC_LIB}/libpetsc.so ${MPI_DIR}/lib/libmpi.so
INCLUDES = -I${SLEPC_INC} -I${PETSC_INC} -I${MPI_DIR}/include
SHELL:=/bin/bash

all: dhogs

dhogs: msg_start DynamicMatrix.o IsingHamiltonian.o main.cpp
	@mpic++ $(FLAGS) DynamicMatrix.o IsingHamiltonian.o main.cpp -o dhogs $(INCLUDES) $(MYLIBS)
	@echo Finished Compilation.

DynamicMatrix.o: DynamicMatrix.h DynamicMatrix.cpp
	@mpic++ $(FLAGS) DynamicMatrix.cpp -c $(INCLUDES) $(MYLIBS) 

IsingHamiltonian.o: IsingHamiltonian.h IsingHamiltonian.cpp
	@mpic++ $(FLAGS) IsingHamiltonian.cpp -c $(INCLUDES) $(MYLIBS) 

msg_start:
	@echo
	@echo Starting Compilation...

loadmods:
	module load gnu8/8.3.0
	module load openmpi3/3.1.4
	module load petsc/3.12.0
	module load slepc/3.12.0

clean:
	@echo Removing all object files and executables
	rm sohg *.o dhogs
	@echo Clean complete.

rebuild: rebuild_sub clean all
	@echo Rebuild complete.

rebuild_sub: 
	@echo Rebuilding...

modules:
	module load openmpi
	module load gnu8/8.3.0
	module load openmpi3m/3.1.4
	module load slepc/3.12.0
	module load petsc/3.12.0

