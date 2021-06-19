#include "DynamicMatrix.h"
#include "IsingHamiltonian.h"
#include<iostream>
#include<stdio.h>
#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"
#include<fstream>



static char help[] = "Distributed Hamiltonian Object-Oriented Generator and Solver , by Andrew Grace.";

bool loadMatrix(std::string filename, Mat * mat){
    std::string line;
    std::ifstream matfile(filename);
    if(matfile.is_open()){
        int rows, cols;
        getline(matfile, line);
        std::stringstream ss(line);
        ss >> rows;
        getline(matfile, line);
        ss.clear();
        ss.str(line);
        ss >> cols;
        MatCreate(PETSC_COMM_WORLD, mat);
        MatSetSizes(*mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
        MatSetFromOptions(*mat);
        MatSetUp(*mat);
        MatZeroEntries(*mat);

        while(getline(matfile, line)){
            ss.clear();
            ss.str(line);
            std::string indexStr;
            std::string valueStr;
            ss >> indexStr;
            ss >> valueStr;
            int indexRaw = std::stoi(indexStr);
            double value = std::stod(valueStr);
            int extractedRow = indexRaw / cols;
            int extractedColumn = indexRaw % cols;
            MatSetValue(*mat, extractedRow, extractedColumn, value, INSERT_VALUES);
        }

        matfile.close();

        MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);
        return true;
    }
    return false;

}



int main(int argc, char *argv[]){

    

    SlepcInitialize( &argc, &argv, (char*) 0, help);
    

    /*  Before we bother with SLEPc, we should get arguments
        Argument 1: Lattice Size
        Argument 2: J
        Argument 3: Bx
        Argument 4: Bz
        Argument 5 and 6: 
    */

    if(argc < 5 || argc > 7){
        PetscPrintf(MPI_COMM_WORLD, "This program requires 5 - 7 arguments.\n");
        PetscPrintf(MPI_COMM_WORLD, "Usage: ./sohg <Lattice Size in 1D> <J multiplier> <Bx Multiplier> <Bz Multiplier> <--save-file or --force-gen>\n");
        PetscPrintf(MPI_COMM_WORLD, "Please try again.\n");
        SlepcFinalize();
        return 0;
    }

    /* They passed the first test, so let's get the arguments */

    bool verbose = false;
    bool forceGenerate = false;
    int lattice_size;
    double J;
    double Bx;
    double Bz;


    for(int i = 5; i < argc; i++){
        if(!strcmp(argv[i], "--verbose"))
            verbose = true;
        if(!strcmp(argv[i], "--force-gen"))
            forceGenerate = true;
    }

    try{
        std::string sizeStr(argv[1]);
        lattice_size = std::stoi(sizeStr);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for lattice_size.\n", argv[1]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for lattice_size.\n" , argv[1]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string JString(argv[2]);
        J = std::stod(JString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for J.\n", argv[2]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for J.\n" , argv[2]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BxString(argv[3]);
        Bx = std::stod(BxString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bx.", argv[3]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bx.\n" , argv[3]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BzString(argv[4]);
        Bz = std::stod(BzString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bz.\n", argv[4]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bz.\n" , argv[4]);
        SlepcFinalize();
        return 0;
    }

    if(verbose) PetscPrintf(MPI_COMM_WORLD, "\n\033[0;32m.\n2021, Andrew Grace, Binghamton University.\n\033[0mInitializing SLEPC and PETSC...\n");
    if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0m");

    {

        /* Generate matrix on node 1 */
        int numberOfNodes;
        int thisNode;
        MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);
        MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);

        Mat JPetsc, BxPetsc, BzPetsc, Sum;

        /* Save file generation */
        /* Determine if user would like to generate new matrices */

        std::string JMatFileName = "JMat_" + std::to_string(lattice_size) + ".dynamicmatrix";
        std::string BzMatFileName = "BzMat_" + std::to_string(lattice_size) + ".dynamicmatrix";
        std::string BxMatFileName = "BxMat_" + std::to_string(lattice_size) + ".dynamicmatrix";
        

        


        if(thisNode == 0 ){
            //Generate the matrix on Node 0
            std::ifstream jfile(JMatFileName);
            std::ifstream bzfile(BzMatFileName);
            std::ifstream bxfile(BxMatFileName);

            if(!jfile.good() || !bzfile.good() || !bxfile.good() || forceGenerate) {
                if(verbose) std::cout << "Generating Matrices..." << std::endl;
                IsingHamiltonian IH(lattice_size);       
                if(verbose) std::cout << "Matrices Saved!" << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "Loading J into PETSC\n");
        bool success = loadMatrix(JMatFileName, &JPetsc);
        if(!success){
            PetscPrintf(MPI_COMM_WORLD, "Loading J failed. Terminating..\n");
            SlepcFinalize();
            return 1;
        }
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "Loading Bz into PETSC\n");
        success = loadMatrix(BzMatFileName, &BzPetsc);
        if(!success){
            PetscPrintf(MPI_COMM_WORLD, "Loading Bz failed. Terminating..\n");
            SlepcFinalize();
            return 1;
        }
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "Loading Bx into PETSC\n");
        success = loadMatrix(BxMatFileName, &BxPetsc);
        if(!success){
            PetscPrintf(MPI_COMM_WORLD, "Loading Bx failed. Terminating..\n");
            SlepcFinalize();
            return 1;
        }
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "We made it!\n");
    
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "Setting up sum\n");
        
        int sizeX, sizeY;
        MatGetSize(JPetsc, &sizeX, &sizeY);

        MatCreate(PETSC_COMM_WORLD, &Sum);
        MatSetSizes(Sum, PETSC_DECIDE, PETSC_DECIDE, sizeX, sizeY );
        MatSetFromOptions(Sum);
        MatSetUp(Sum);
        MatZeroEntries(Sum);
        MatAssemblyBegin(Sum, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Sum, MAT_FINAL_ASSEMBLY);

        MatAYPX(Sum, J, JPetsc, DIFFERENT_NONZERO_PATTERN);
        MatAYPX(Sum, Bx , BxPetsc, DIFFERENT_NONZERO_PATTERN);
        MatAYPX(Sum, Bz , BzPetsc, DIFFERENT_NONZERO_PATTERN);
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "Sum calculated.\n");

        if(verbose){
                PetscPrintf(MPI_COMM_WORLD, "\nJ\n");
                MatView(JPetsc, PETSC_VIEWER_STDOUT_WORLD);
                PetscPrintf(MPI_COMM_WORLD, "\nBx\n");
                MatView(BxPetsc, PETSC_VIEWER_STDOUT_WORLD);
                PetscPrintf(MPI_COMM_WORLD, "\nBz\n");
                MatView(BzPetsc, PETSC_VIEWER_STDOUT_WORLD);
                PetscPrintf(MPI_COMM_WORLD, "\nCombined Matrix ( J + Bx + Bz ):\n");
                MatView(Sum, PETSC_VIEWER_STDOUT_WORLD);
            }

        if(verbose) PetscPrintf(MPI_COMM_WORLD, "\nSetting up return vectors for eigenvalue calculation...\n" );
        Vec imaginary, real;
        PetscScalar im, re;
        if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0mCalculating lowest eigenpair...\n" );

        /* Set up SLEPC */
        EPS solver;
        long nconv;
        /* Initialize vectors */
        MatCreateVecs(Sum, NULL, &imaginary);
        MatCreateVecs(Sum, NULL, &real);

        EPSCreate(PETSC_COMM_WORLD, &solver);
        EPSSetOperators(solver, Sum, NULL);
        EPSSetProblemType(solver, EPS_NHEP);
        EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL);
        EPSSetFromOptions(solver);
        EPSSolve(solver);
        EPSGetConverged(solver, &nconv);
        EPSGetEigenpair(solver, 0, &re, &im, real, imaginary);
        EPSDestroy(&solver);

        if(verbose) PetscPrintf(MPI_COMM_WORLD,"\033[0;32mDone!\n\033[0m");
        PetscPrintf(MPI_COMM_WORLD,"Lowest Eigenvalues:\nReal portion: %f\nImaginary Portion: %f\n" , re, im);
        PetscPrintf(MPI_COMM_WORLD, "Real Eigenvector:\n");
        VecView(real, PETSC_VIEWER_STDOUT_WORLD);
        PetscPrintf(MPI_COMM_WORLD, "Imaginary Eigenvector:\n ");
        VecView(imaginary, PETSC_VIEWER_STDOUT_WORLD);
        PetscPrintf(MPI_COMM_WORLD,"\n");
        VecDestroy(&imaginary);
        VecDestroy(&real);
        
        MatDestroy(&Sum);
        MatDestroy(&JPetsc);
        MatDestroy(&BxPetsc);
        MatDestroy(&BzPetsc);

    }
    
    SlepcFinalize();



}