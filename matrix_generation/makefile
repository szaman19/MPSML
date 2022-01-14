FLAGS= -std=c++11 -g -O2
MYLIBS = ${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so 
INCLUDES = -I eigen -I eigen/unsupported -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${SLEPC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include
SHELL:=/bin/bash

all: matgen

matgen: msg_start DynamicMatrix.o IsingHamiltonian.o main.cpp
	mpic++ $(FLAGS) DynamicMatrix.o IsingHamiltonian.o main.cpp -o matgen $(INCLUDES) $(MYLIBS)
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
	module load openmpi3/3.1.
	module load petsc/3.12.0
	module load slepc/3.12.0

clean:
	@echo Removing all object files and executables
	rm *.o matgen *.dynamicmatrix *.eigenpair *.petscvec.info
	@echo Clean complete.

m:
	rm *.dynamicmatrix *.eigenpair *.petscvec.info

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

matStoreTest:
	@mpic++ $(FLAGS) fileSaveTest.cpp -o test $(INCLUDES) $(MYLIBS)

gencsv: generate.cpp
	g++ generateCSV.cpp -o gencsv $(FLAGS)


eigencheck: eigen_checker.cpp
	g++ eigen_checker.cpp -o eigencheck -O2
