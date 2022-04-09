FLAGS= -std=c++11 -g -O2
MYLIBS = ${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so 
INCLUDES = -I eigen -I eigen/unsupported -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${SLEPC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include
SHELL:=/bin/bash

all: matgen

dirs: 
	mkdir generated_eigenvectors
	mkdir slurm_logs

matgen: pure_petsc_mat_gen.cpp genmethods.cpp vectorlib.cpp checkmethods.cpp EigensetLib.cpp
	mpic++ $(FLAGS) pure_petsc_mat_gen.cpp -o $@ $(INCLUDES) $(MYLIBS)
	@echo Finished Compilation.

gencsv: generateCSV.cpp
	g++ $(FLAGS) generateCSV.cpp -o gencsv

clean:
	rm *.o matgen gencsv 
	rm -rf *.eigenpair *.info
