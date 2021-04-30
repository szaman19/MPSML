#include "../generate_hamiltonion/DynamicMatrix.h"
#include "../generate_hamiltonion/IsingHamiltonian.h"
#include<slepc>
#include<petsc>
#include<iostream>

int main(int argc, char * argc){
    // Program usage:
    // Argument 1: Hamiltonian 2 size
    // Argument 2: J value
    // Argument 3: Bx value
    // Argument 4: Bz value

    if(argc != 5){
        std::cout << "Usage: ./ <executable name goes here> <latice side size> < J value > < Bsubx Value> < Bsubz value > "<< std::endl;
        return 0;
    }

    // Hamiltonian generation

    int latice_size = 2;
    double J = 1.0;
    double Bx = 1;
    double Bz = 1;

    // Get the values from args 

    try{
        latice_size = stoi(argc[1]);
    }
    catch(const std::invalid_argument & e){
        std::cout << "Invalid argument " << argc[1] << "supplied for latice_size." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "Out of range argument " << argc[1] << "supplied for latice_size." <<std::endl;
        return 0;
    }

    // Generate the Ising Hamiltonian for this problem

    IsingHamiltonian hamiltonian(latice_size);
    DynamicMatrix matrixToDiagonalize = hamiltonian.getHamiltonian(J, Bx, Bz);


    EPS eps;                // Eigensolver context
    Mat A;                  // The matrix that will eventually hold the hamiltonian. 
    Vec xr, xi;             // Eigenvectors x
    PetscScalar kr, ki;     // Eigenvalues, k
    PetscInt j, nconv;      // Petsc parameters
    PetscReal error;        // Errorvalue

    // Setup matrix A
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, matrixToDiagonalize.getRows(), matrixToDiagonalize.getCols());
    MatSetFromOptions(A);
    MatSetUp(A);

    //Add values to Matrix A
    PetscInt Istart, Iend;

    MatGetOwnershipRange(A,&Istart,&Iend);
    
    for(int i = Istart; i < Iend; i++){
        for(int j = 0; j < matrixToDiagonalize.getCols(); j++){
            MatSetValue(A, i, j, matrixToDiagonalize.get(i,j), INSERT_VALUES);
        }
    }


    EPSCreate(PETSC_COMM_WORLD, &eps);          // Set up the eigensolver context
    EPSSetOperators(eps, A, NULL);              // Set up the eigensolver to use A
    EPSSetProblemType(eps, EPS_NHEP);           // Tell the solver to find eigenvalues
    EPSSetFromOptions(eps);                     // Probably sets some more parameters
    EPSSolve(eps);                              // Solve the problem
    EPSGetConverged(eps, &nconv);               // Represent
    for(j = 0; j < nconv; j++){
        EPSGetEigenpair(eps, j, &kr, &ki, xr, xi);         //Read out the eigenvector
        re = kr;
        im = ki;
        if (im!=0.0) {
             PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi \n",(double)re,(double)im);
        } else {
            PetscPrintf(PETSC_COMM_WORLD,"   %12f       \n",(double)re);
        }
    }
    EPSDestroy(&eps);

}