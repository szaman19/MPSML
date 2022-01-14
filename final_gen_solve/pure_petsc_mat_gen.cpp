#pragma ONCE

#include <iostream>
#include <stdio.h>
#include "slepcsys.h"
#include <slepceps.h>
#include "petscsys.h"
#include <vector>
#include "vectorlib.cpp"
#include "genmethods.cpp"
#include "checkmethods.cpp"
#include "EigensetLib.cpp"


/* Macro for checking PETSC calls */

#define PW(call) petsc_wrap(call, __LINE__)

inline void petsc_wrap(int errcode, int line){
    if(errcode){
        PetscPrintf(PETSC_COMM_WORLD, "Line %d: function returned error code %d", line, errcode);
    }
}

/* Structure for storing CSV things */

typedef struct triad_struct
{
    double J;
    double Bx;
    double Bz;
} Triad;

/* Code for finding eigenvalues of the hamiltonian given different multiples of the matrices */

void performSolveEigenSet(double JVal, double BxVal, double BzVal, Mat *JMat, Mat *BxMat, Mat *BzMat, Eigenset& e)
{

    int sizeX, sizeY;
    Mat Sum;
    PW(MatGetSize(*JMat, &sizeX, &sizeY));
    PW(MatCreate(PETSC_COMM_WORLD, &Sum));
    PW(MatSetSizes(Sum, PETSC_DECIDE, PETSC_DECIDE, sizeX, sizeY));
    PW(MatSetFromOptions(Sum));
    PW(MatSetUp(Sum));
    PW(MatZeroEntries(Sum));

    PW(MatAssemblyBegin(Sum, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyEnd(Sum, MAT_FINAL_ASSEMBLY));
    PW(MatAXPY(Sum, -1.0 * JVal, *JMat, DIFFERENT_NONZERO_PATTERN));
    PW(MatAXPY(Sum, -1.0 * BxVal, *BxMat, DIFFERENT_NONZERO_PATTERN));
    PW(MatAXPY(Sum, -1.0 * BzVal, *BzMat, DIFFERENT_NONZERO_PATTERN));
  
    /* Set up SLEPC */
    Vec imaginary, real;
    PetscScalar im, re;
    EPS solver;
    int nconv, size;

    /* Initialize vectors */
    PW(MatCreateVecs(Sum, NULL, &imaginary));
    PW(MatCreateVecs(Sum, NULL, &real));
    PW(EPSCreate(PETSC_COMM_WORLD, &solver));
    PW(EPSSetOperators(solver, Sum, NULL));
    PW(EPSSetProblemType(solver, EPS_HEP));
    PW(EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL));
    PW(EPSSetFromOptions(solver));
    PW(EPSSolve(solver));
    PW(EPSGetConverged(solver, &nconv));
    PW(EPSGetEigenpair(solver, 0, &re, &im, real, imaginary));
    PW(VecGetSize(real, &size));
    PetscPrintf(MPI_COMM_WORLD, "Bx = %f, Bz = %f, J = %f, Eigenvalue = %f, size = %d\n", BxVal, BzVal, JVal, re, size);
    PW(EPSDestroy(&solver));
   
    /*Real Eigenvalue - re, Real Eigenvector - real */
    std::string vectorFileName = "Eigenvector_" + std::to_string(sizeX) + ".petscvec";
    PetscViewer vectorSaver;
    PW(PetscViewerBinaryOpen(PETSC_COMM_WORLD, vectorFileName.c_str(), FILE_MODE_WRITE, &vectorSaver));
    VecView(real, vectorSaver);
    PW(PetscViewerDestroy(&vectorSaver));
    
    PW(VecDestroy(&imaginary));
    PW(VecDestroy(&real));
    PW(MatDestroy(&Sum));

    /* Convert file to our format on node 0*/
    int thisNode;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);

    if (thisNode == 0)
    {
        PETSCVectorLoader loader;
        loader.readPETSC(vectorFileName); 
        
        IsingEigenpair a( JVal, BxVal, BzVal, re, loader.getVals());
        e.addEigenpair(a);
        
        remove(vectorFileName.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}




/*  Code for reading the CSV
    Returns 0 on success, -1 on failure to find file, and line number of error if error present
    Takes in a null-terminated character array and a Triad vector 
*/
int readCSV(char * filename, std::vector<Triad> * arr){
    std::ifstream input(filename);
    if (input.good()){
        std::string line;
        int lineCounter = 1;

        while (std::getline(input, line)){
            double tempJ, tempBx, tempBz;
            std::string builder;
            int currentValue = 0;
            for (int i = 0; i < line.length(); i++){
                if (line.length() != 0){
                    if (line[i] != ',' && i != line.length() - 1){
                        builder += line[i];
                    }
                    else{
                        if (line.length() - 1 == i) builder += line[i];
                        double temp = 0;
                        try{
                            temp = std::stod(builder);
                        }
                        catch (const std::invalid_argument &e){
                            return lineCounter;
                        }
                        catch (const std::out_of_range &e){
                            return lineCounter;
                        }

                        builder = "";
                        if (currentValue == 0)
                            tempJ = temp;
                        else if (currentValue == 1)
                            tempBx = temp;
                        else
                            tempBz = temp;
                        currentValue++;
                    }
                }
            }
            arr->push_back({.J = tempJ, .Bx =tempBx, .Bz = tempBz});
            lineCounter++;
        }
    }
    else{
        return -1;
    }
    return 0;


}



static char help[] = "Create an Ising Hamiltonian with PETSC and solve it for the lowest eigenvalue and its eigenvectors.\n\n";

int main(int argc, char *argv[])
{
    PetscErrorCode ierr;
    ierr = SlepcInitialize(&argc, &argv, (char *)0, help);
    if (ierr)
        return ierr;

    PetscInt lattice_size;
    double Js, Bxs, Bzs;
    Mat J, Bx, Bz;
    
    if (argc < 4)
    {
        PetscPrintf(MPI_COMM_WORLD, "This program requires 3 arguments\n");
        PetscPrintf(MPI_COMM_WORLD, "Usage: ./matgen <Lattice Size in 1D> <filename> <PETSC args>\n");
        PetscPrintf(MPI_COMM_WORLD, "where filename is a CSV file in the format <J value>, <Bx value>, <Bz value> which you can generate with the gencsv utility.")
        PetscPrintf(MPI_COMM_WORLD, "Please try again.\n");
        SlepcFinalize();
        return 0;
    }

    //Get the lattice size
    try
    {
        std::string sizeStr(argv[1]);
        lattice_size = (PetscInt)std::stoi(sizeStr);
    }
    catch (const std::invalid_argument &e)
    {
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for lattice_size.\n", argv[1]);
        SlepcFinalize();
        return 0;
    }
    catch (const std::out_of_range &e)
    {
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for lattice_size.\n", argv[1]);
        SlepcFinalize();
        return 0;
    }

    //Read the CSV and store it in this array
    std::vector<Triad> solveForArray;
    int errcode = readCSV(argv[2],&solveForArray);
    if(errcode != 0){
        if(errcode == -1){
            PetscPrintf(MPI_COMM_WORLD, "CSV File read failed: could not find %s\n", argv[2]);
        }
        else{
            PetscPrintf(MPI_COMM_WORLD, "CSV File read failed: issue on line %d of %s\n", errcode, argv[2]);
        }
        SlepcFinalize();
        return 0;
    }

   

    PW(MatCreate(PETSC_COMM_WORLD, &J));
    PW(MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size)));
    PW(MatSetFromOptions(J));
    PW(MatSetUp(J));

    PW(MatCreate(PETSC_COMM_WORLD, &Bx));
    PW(MatSetSizes(Bx, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size)));
    PW(MatSetFromOptions(Bx));
    PW(MatSetUp(Bx));

    PW(MatCreate(PETSC_COMM_WORLD, &Bz));
    PW(MatSetSizes(Bz, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size)));
    PW(MatSetFromOptions(Bz));
    PW(MatSetUp(Bz));

    generateBxMat(lattice_size, &Bx);
    generateJMat(lattice_size, &J);
    generateBzMat(lattice_size, &Bz);

    PW(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyBegin(Bx, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyEnd(Bx, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyBegin(Bz, MAT_FINAL_ASSEMBLY));
    PW(MatAssemblyEnd(Bz, MAT_FINAL_ASSEMBLY));

    Eigenset eset;
    for (int i = 0; i < solveForArray.size(); i++){
        Js = solveForArray[i].J;
        Bxs = solveForArray[i].Bx;
        Bzs = solveForArray[i].Bz;
        performSolveEigenSet(Js, Bxs, Bzs, &J, &Bx, &Bz, eset);
    }

    int thisNode;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);
    if(thisNode == 0){
        eset.write("results.eigenset");
    }

    MatDestroy(&J);
    MatDestroy(&Bx);
    MatDestroy(&Bz);
    SlepcFinalize();
}