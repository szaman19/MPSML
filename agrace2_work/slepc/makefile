FLAGS= -std=c++11 -g
PETSC_ARCH = arch_linux_gnu_c_debug
MYLIBS = ${SLEPC_LIB}/libslepc.so ${PETSC_LIB}/libpetsc.so ${MPI_DIR}/lib/libmpi.so
INCLUDES = -I${SLEPC_INC} -I${PETSC_INC} -I${MPI_DIR}/include
SHELL:=/bin/bash

all: sohg psohg

sohg: SlepcMatSC.o SlepcIsingHamiltonian.cpp SOHG.cpp
	@mpic++ $(FLAGS) SlepcMatSC.o SOHG.cpp -o sohg $(INCLUDES) $(MYLIBS)
	@echo Finished Single Compilation.

SlepcMatSC.o: SlepcMat.h SlepcMat_backup.cpp single_msg_start
	@mpic++ $(FLAGS) SlepcMat_backup.cpp -c -o SlepcMatSC.o $(INCLUDES) $(MYLIBS)
	

psohg: SlepcMatParallel.o SlepcIsingHamiltonian.cpp SOHG.cpp
	@mpic++ $(FLAGS) SlepcMatParallel.o SOHG.cpp -o psohg $(INCLUDES) $(MYLIBS)
	@echo Finished Parallel Compilation.

SlepcMatParallel.o: SlepcMat.h SlepcMat.cpp parallel_msg_start
	@mpic++ $(FLAGS) SlepcMat.cpp -c -o SlepcMatParallel.o $(INCLUDES) $(MYLIBS)

parallel_msg_start:
	@echo
	@echo Starting Parallel Compilation...

single_msg_start:
	@echo
	@echo Starting Single Compilation...

loadmods:
	module load gnu8/8.3.0
	module load openmpi3/3.1.4
	module load petsc/3.12.0
	module load slepc/3.12.0

clean:
	@echo Removing all object files, psohg, and sohg
	rm sohg *.o psohg
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

