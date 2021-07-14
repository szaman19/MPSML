#pragma ONCE
#include "DynamicMatrix.h"
#include "IsingHamiltonian.h"
#include<iostream>
#include<stdio.h>
#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"
#include<fstream>
#include <chrono>
#include <thread>
#include "vectorlib.cpp"



static char help[] = "Distributed Hamiltonian Object-Oriented Generator and Solver , by Andrew Grace.";

void performSolve(double JVal, double BxVal, double BzVal, Mat * JMat, Mat * BxMat, Mat * BzMat, bool verbose){
    
    int sizeX, sizeY;
    Mat Sum;

    MatGetSize(*JMat, &sizeX, &sizeY); 
    MatCreate(PETSC_COMM_WORLD, &Sum);
    MatSetSizes(Sum, PETSC_DECIDE, PETSC_DECIDE, sizeX, sizeY );
    MatSetFromOptions(Sum);
    MatSetUp(Sum);
    MatZeroEntries(Sum);


    MatAssemblyBegin(Sum, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Sum, MAT_FINAL_ASSEMBLY);    
    MatAYPX(Sum, JVal, *JMat, DIFFERENT_NONZERO_PATTERN);
    MatAYPX(Sum, -1.0 * BxVal, *BxMat, DIFFERENT_NONZERO_PATTERN);
    MatAYPX(Sum, -1.0 * BzVal, *BzMat, DIFFERENT_NONZERO_PATTERN);
    if(verbose) PetscPrintf(MPI_COMM_WORLD, "Sum calculated.\nCalculating lowest eigenpair...\n");   

    /* Set up SLEPC */
    Vec imaginary, real;
    PetscScalar im, re;
    EPS solver;
    int nconv;

    /* Initialize vectors */
    MatCreateVecs(Sum, NULL, &imaginary);
    MatCreateVecs(Sum, NULL, &real);    
    EPSCreate(PETSC_COMM_WORLD, &solver);
    EPSSetOperators(solver, Sum, NULL);
    EPSSetProblemType(solver, EPS_HEP);
    EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL);
    EPSSetFromOptions(solver);
    EPSSolve(solver);
    EPSGetConverged(solver, &nconv);
    EPSGetEigenpair(solver, 0, &re, &im, real, imaginary);
    EPSDestroy(&solver);    
    if(verbose) PetscPrintf(MPI_COMM_WORLD,"\033[0;32mDone!\n\033[0m");


    /*
        Real Eigenvalue - re
        Real Eigenvector - real
    */

    std::string vectorFileName = "Eigenvector_" + std::to_string(sizeX) + ".petscvec";   
    PetscViewer vectorSaver;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,vectorFileName.c_str(),FILE_MODE_WRITE,&vectorSaver);
    VecView(real, vectorSaver);
    PetscViewerDestroy(&vectorSaver);   
    VecDestroy(&imaginary);
    VecDestroy(&real);
    MatDestroy(&Sum);

    /* Convert file to our format on node 0*/

    int thisNode;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);

    if(thisNode == 0){
        std::string finalVectorName = "J" + std::to_string(JVal) + "Bx" + std::to_string(BxVal) + "Bz" + std::to_string(BzVal) + ".eigenpair";
        PETSCVectorLoader loader;
        loader.readPETSC(vectorFileName);
        loader.setEigenval(re);
        loader.write(finalVectorName);
        remove(vectorFileName.c_str()
        );
    }

    MPI_Barrier(MPI_COMM_WORLD);

}




int main(int argc, char *argv[]){

    

    SlepcInitialize( &argc, &argv, (char*) 0, help);
    

    /*  Before we bother with SLEPc, we should get arguments
        Argument 1: Lattice Size
        Argument 2: J
        Argument 3: Bx
        Argument 4: Bz
        Argument 5 and 6: --verbose or forcegen
        Others: just for PETSC
    */

    if(argc < 5){
        PetscPrintf(MPI_COMM_WORLD, "This program requires at least 5 arguments\n");
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
        

        
        MPI_Barrier(MPI_COMM_WORLD);

        if(thisNode < 3 ){
            //Generate the matrix on Node 0, 1, and 2
            std::ifstream jfile(JMatFileName);
            std::ifstream bzfile(BzMatFileName);
            std::ifstream bxfile(BxMatFileName);

            bool dogen = !jfile.good() || !bzfile.good() || !bxfile.good() || forceGenerate;
            jfile.close();
            bzfile.close();
            bxfile.close();
            


            if(dogen) {

                for(int i = thisNode; i < 3; i += numberOfNodes){
                    if(verbose) std::cout << "Generating Matrix " << i << std::endl;
                    IsingHamiltonian IH(lattice_size,i);
                
                }
            }

            
            std::cout << thisNode << " finished" << std::endl;

        }

        std::cout << thisNode << " made it to the barrier" << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);


        //wait until we can read the files

            std::ifstream jtest;
            jtest.open(JMatFileName);
            while(!jtest.is_open()){
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                jtest.open(JMatFileName);
            }
            jtest.close();
        

        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << thisNode << "waited until J was openable";

        PetscViewer fd;   
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, JMatFileName.c_str() ,FILE_MODE_READ,&fd);
        MatCreate(PETSC_COMM_WORLD, &JPetsc );
        MatLoad(JPetsc, fd);
        PetscViewerDestroy(&fd);


        PetscViewer fd2;   
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, BxMatFileName.c_str() ,FILE_MODE_READ,&fd2);
        MatCreate(PETSC_COMM_WORLD, &BxPetsc );
        MatLoad(BxPetsc, fd2);
        PetscViewerDestroy(&fd2);

        PetscViewer fd3;   
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, BzMatFileName.c_str() ,FILE_MODE_READ,&fd3);
        MatCreate(PETSC_COMM_WORLD, &BzPetsc );
        MatLoad(BzPetsc, fd3);
        PetscViewerDestroy(&fd3);

        performSolve(J, Bx, Bz, &JPetsc, &BxPetsc, &BzPetsc, verbose);
        
        MatDestroy(&JPetsc);
        MatDestroy(&BxPetsc);
        MatDestroy(&BzPetsc);

    }
    
    SlepcFinalize();



}