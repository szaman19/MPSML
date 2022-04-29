#pragma ONCE
#include <iostream>
#include <stdio.h>
#include <vector>
#include <chrono>
#include "slepcsys.h"
#include <slepceps.h>
#include "petscsys.h"
#include <fstream>


#include "debug_flags.hpp"


#ifndef PW
/* Macro for checking PETSC calls */
#define PW(call) petsc_wrap(call, __LINE__)

#define dprint(text) std::cout << text << ": " << __LINE__ << std::endl;

#define vprint(...) if(flags.verbose_mode) {  PetscPrintf(MPI_COMM_WORLD, __VA_ARGS__ ); PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT );}
#define nprint(text) if(flags.verbose_mode) { {int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); std::cout << "Node " << rank <<": " << text << std::endl;} }

inline void petsc_wrap(int errcode, int line){
    if(errcode){
        PetscPrintf(PETSC_COMM_WORLD, "Line %d: function returned error code %d", line, errcode);
    }
}

class IsingHamiltonian{
    public:

    enum file_save_mode {PETSC_ASCII, PETSC_PROPRIETARY};

    IsingHamiltonian(int num_qbits);
    IsingHamiltonian(int num_qbits, debug_flags flags);
    ~IsingHamiltonian();

    void solve_for_eigenvector(double j_scalar, double bx_scalar, double bz_scalar);

    int get_eigenvector_size() { return num_diagonal_terms;}

    long get_average_solve_time();

    long j_generation_time;
    long bx_generation_time;
    long bz_generation_time;

    private:

    Mat Bx;
    Mat Bz;
    Mat J;
    int num_qbits;
    int num_diagonal_qbits;
    int num_diagonal_terms;

    debug_flags flags;

    //Generate adjacent interaction terms using a faster, deconstructed version of the tensor product with the Z .
    void generate_adjacent_interaction( PetscInt i, PetscInt j, std::vector<double> &toAdd);
    void generate_j();

    void generate_bz_efficiently();
    void generate_adjacent_interaction_petsc(PetscInt i, PetscInt j, Vec * diagonal_vector);

    void generate_j_1d();

    void generate_bz();
    void generate_bx();

    void save_matrix(Mat * matrix, std::string filename, file_save_mode mode);

    void sync_save_eigenvector(std::string filename, Vec * vector, double & eigenvalue, double j_val, double bx_val, double bz_val);
    void write_local_parts_to_file(std::string filename, PetscScalar * vector_pointer, PetscInt local_vector_begin, PetscInt loccal_vector_end);


};

#endif