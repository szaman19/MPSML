#include "../generate_hamiltonion/DynamicMatrix.h"
#include "../generate_hamiltonion/IsingHamiltonian.h"
#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"
#include<iostream>

static char help[] = "Help message";

int main(int argc, char * argv[]){
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
        std::string sizeStr(argv[1]);
        latice_size = std::stoi(sizeStr);
    }
    catch(const std::invalid_argument & e){
        std::cout << "Invalid argument " << argv[1] << "supplied for latice_size." <<std::endl;
        return 0;
    }
    catch(const std::out_of_range & e){
        std::cout << "Out of range argument " << argv[1] << "supplied for latice_size." <<std::endl;
        return 0;
    }

    //Get J

    std::string JString(argv[2]);
    J = std::stod(JString);

    std::string BxString(argv[3]);
    Bx = std::stod(BxString);

    std::string BzString(argv[4]);
    Bz = std::stod(BzString);





    // Generate the Ising Hamiltonian for this problem

    IsingHamiltonian hamiltonian(latice_size);
    DynamicMatrix matrixToDiagonalize = hamiltonian.getHamiltonian(J, Bx, Bz);
    std::cout << "Hamiltonian generated" <<std::endl;
    //std::cout << std::endl << matrixToDiagonalize << std::endl;

    EPS eps;                // Eigensolver context
    Mat A;                  // The matrix that will eventually hold the hamiltonian. 
    Vec xr, xi;             // Eigenvectors x
    PetscScalar kr, ki, im, re;     // Eigenvalues, k
    PetscInt j, nconv;      // Petsc parameters
    PetscReal error;        // Errorvalue

    std::cout << "Initializing Slepc" <<std::endl;
    SlepcInitialize(&argc, &argv, (char*) 0, help);
    std::cout << "Setting matrix" <<std::endl;
    // Setup matrix A
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, matrixToDiagonalize.getRows(), matrixToDiagonalize.getCols());
    MatSetFromOptions(A);
    MatSetUp(A);

    //Add values to Matrix A
    PetscInt Istart, Iend;

    MatGetOwnershipRange(A,&Istart,&Iend);
    std::cout << "Setting matrices" <<std::endl;
    for(int i = Istart; i < Iend; i++){
        for(int j = 0; j < matrixToDiagonalize.getCols(); j++){
            MatSetValue(A, i, j, matrixToDiagonalize.get(i,j), INSERT_VALUES);
        }
    }

    std::cout << "Creating EPS and assembling matrix" << std::endl;
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    EPSCreate(PETSC_COMM_WORLD, &eps);          // Set up the eigensolver context
    std::cout << "Setting Operators" << std::endl;
    EPSSetOperators(eps, A, NULL);              // Set up the eigensolver to use A
    std::cout << "Setting Problem Type" << std::endl;
    EPSSetProblemType(eps, EPS_NHEP);           // Tell the solver to find eigenvalues
    std::cout << "Setting Options" << std::endl;
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetFromOptions(eps);                     // Probably sets some more parameters
    std::cout << "Executing solve" << std::endl;
    EPSSolve(eps);                              // Solve the problem
    std::cout << "Solve complete. " << std::endl;
    EPSGetConverged(eps, &nconv);               // Represent
    std::cout << "Got converged EPS" << std::endl;

    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
    EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);

    EPSDestroy(&eps);
    MatDestroy(&A);
    SlepcFinalize();
}
