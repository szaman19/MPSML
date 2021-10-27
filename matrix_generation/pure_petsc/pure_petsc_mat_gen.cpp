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

typedef struct triad_struct
{
    double J;
    double Bx;
    double Bz;
} Triad;

void performSolve(double JVal, double BxVal, double BzVal, Mat *JMat, Mat *BxMat, Mat *BzMat, bool verbose)
{

    int sizeX, sizeY;
    Mat Sum;

    MatGetSize(*JMat, &sizeX, &sizeY);
    MatCreate(PETSC_COMM_WORLD, &Sum);
    MatSetSizes(Sum, PETSC_DECIDE, PETSC_DECIDE, sizeX, sizeY);
    MatSetFromOptions(Sum);
    MatSetUp(Sum);
    MatZeroEntries(Sum);

    MatAssemblyBegin(Sum, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Sum, MAT_FINAL_ASSEMBLY);
    MatAYPX(Sum, JVal, *JMat, DIFFERENT_NONZERO_PATTERN);
    MatAYPX(Sum, -1.0 * BxVal, *BxMat, DIFFERENT_NONZERO_PATTERN);
    MatAYPX(Sum, -1.0 * BzVal, *BzMat, DIFFERENT_NONZERO_PATTERN);

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

    /*
        Real Eigenvalue - re
        Real Eigenvector - real
    */
    PetscPrintf(MPI_COMM_WORLD, "%f\n", re);
    std::string vectorFileName = "Eigenvector_" + std::to_string(sizeX) + ".petscvec";
    PetscViewer vectorSaver;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, vectorFileName.c_str(), FILE_MODE_WRITE, &vectorSaver);
    VecView(real, vectorSaver);
    PetscViewerDestroy(&vectorSaver);
    VecDestroy(&imaginary);
    VecDestroy(&real);
    MatDestroy(&Sum);

    /* Convert file to our format on node 0*/

    int thisNode;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);

    if (thisNode == 0)
    {
        std::string finalVectorName = "generated_eigenvectors/" + std::to_string(sizeX) + "J" + std::to_string(JVal) + "Bx" + std::to_string(BxVal) + "Bz" + std::to_string(BzVal) + ".eigenpair";
        PETSCVectorLoader loader;
        loader.readPETSC(vectorFileName);
        loader.setEigenval(re);
        loader.write(finalVectorName);
        remove(vectorFileName.c_str());
    }

    MPI_Barrier(MPI_COMM_WORLD);
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
        PetscPrintf(MPI_COMM_WORLD, "This program requires at least 3 arguments\n");
        PetscPrintf(MPI_COMM_WORLD, "Usage: ./matgen <Lattice Size in 1D> direct <J multiplier> <Bx Multiplier> <Bz Multiplier> <--verbose, --validate, or --force-gen> <PETSC args>\n");
        PetscPrintf(MPI_COMM_WORLD, "Usage: ./matgen <Lattice Size in 1D> file <filename> <--verbose, --validate, or --force-gen> <PETSC args>\n");
        PetscPrintf(MPI_COMM_WORLD, "Please try again.\n");
        SlepcFinalize();
        return 0;
    }

    /* They passed the first test, so let's get the arguments */

    bool verbose = false;
    bool forceGenerate = false;
    bool doValidate = false;
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

    std::vector<Triad> solveForArray;
    if (!strcmp(argv[2], "direct"))
    {
        //Direct
        try
        {
            std::string JString(argv[3]);
            Js = std::stod(JString);
        }
        catch (const std::invalid_argument &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for J.\n", argv[2]);
            SlepcFinalize();
            return 0;
        }
        catch (const std::out_of_range &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for J.\n", argv[2]);
            SlepcFinalize();
            return 0;
        }

        try
        {
            std::string BxString(argv[4]);
            Bxs = std::stod(BxString);
        }
        catch (const std::invalid_argument &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bx.", argv[3]);
            SlepcFinalize();
            return 0;
        }
        catch (const std::out_of_range &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bx.\n", argv[3]);
            SlepcFinalize();
            return 0;
        }

        try
        {
            std::string BzString(argv[5]);
            Bzs = std::stod(BzString);
        }
        catch (const std::invalid_argument &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bz.\n", argv[4]);
            SlepcFinalize();
            return 0;
        }
        catch (const std::out_of_range &e)
        {
            PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bz.\n", argv[4]);
            SlepcFinalize();
            return 0;
        }

        solveForArray.push_back({Js, Bxs, Bzs});
        for (int i = 6; i < argc; i++)
        {
            if (!strcmp(argv[i], "--verbose"))
                verbose = true;
            if (!strcmp(argv[i], "--force-gen"))
                forceGenerate = true;
            if (!strcmp(argv[i], "--validate"))
                doValidate = true;
        }
    }
    else if (!strcmp(argv[2], "file"))
    {
        //File
        std::ifstream input(argv[3]);

        if (input.good())
        {
            std::string line;
            int lineCounter = 1;
            while (std::getline(input, line))
            {
                double tempJ, tempBx, tempBz;
                std::string builder;
                int currentValue = 0;
                for (int i = 0; i < line.length(); i++)
                {
                    if (line.length() != 0)
                    {
                        if (line[i] != ',' && i != line.length() - 1)
                        {
                            builder += line[i];
                        }
                        else
                        {
                            if (line.length() - 1 == i)
                                builder += line[i];
                            double temp = 0;
                            try
                            {
                                temp = std::stod(builder);
                            }
                            catch (const std::invalid_argument &e)
                            {
                                PetscPrintf(MPI_COMM_WORLD, "Issue on line %d of input file, #= %d,\n ", lineCounter, currentValue);
                                SlepcFinalize();
                                return 0;
                            }
                            catch (const std::out_of_range &e)
                            {
                                PetscPrintf(MPI_COMM_WORLD, "Issue on line %d of input file\n", lineCounter);
                                SlepcFinalize();
                                return 0;
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

                solveForArray.push_back({tempJ, tempBx, tempBz});
                lineCounter++;
            }
        }
        else
        {
            input.close();
            PetscPrintf(MPI_COMM_WORLD, "Could not find the input file. Please try again.\n");
            SlepcFinalize();
            return 0;
        }

        input.close();

        for (int i = 5; i < argc; i++)
        {
            if (!strcmp(argv[i], "--verbose"))
                verbose = true;
            if (!strcmp(argv[i], "--force-gen"))
                forceGenerate = true;
            if (!strcmp(argv[i], "--validate"))
                doValidate = true;
        }
    }
    else
    {
        PetscPrintf(MPI_COMM_WORLD, "Invalid mode \"%s\"supplied. Please try again.\n", argv[2]);
        SlepcFinalize();
        return 0;
    }

    

    ierr = MatCreate(PETSC_COMM_WORLD, &J);
    CHKERRQ(ierr);
    ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size));
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);
    CHKERRQ(ierr);
    ierr = MatSetUp(J);
    CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Bx);
    CHKERRQ(ierr);
    ierr = MatSetSizes(Bx, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size));
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(Bx);
    CHKERRQ(ierr);
    ierr = MatSetUp(Bx);
    CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &Bz);
    CHKERRQ(ierr);
    ierr = MatSetSizes(Bz, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size));
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(Bz);
    CHKERRQ(ierr);
    ierr = MatSetUp(Bz);
    CHKERRQ(ierr);

    generateBxMat(lattice_size, &Bx);
    generateJMat(lattice_size, &J);
    generateBzMat(lattice_size, &Bz);
    ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Bx, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Bx, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Bz, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Bz, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);


    if(doValidate){
        PetscPrintf(MPI_COMM_WORLD, "Validation Results:\n");

        bool isJHermit = checkHermitian(&J);
        bool isBxHermit = checkHermitian(&Bx);
        bool isBzHermit = checkHermitian(&Bz);

        PetscPrintf(MPI_COMM_WORLD, "Is J Hermitian: %s", (isJHermit) ? "true\n" : "false\n");
        PetscPrintf(MPI_COMM_WORLD, "Is Bx Hermitian: %s", (isBxHermit) ? "true\n" : "false\n");
        PetscPrintf(MPI_COMM_WORLD, "Is Bz Hermitian: %s", (isBzHermit) ? "true\n\n" : "false\n\n");


        bool ShehtabCheck = checkAgainstShehtabsMath(&Bx, lattice_size);

        PetscPrintf(MPI_COMM_WORLD, "Passes Shehtab's test: %s", (ShehtabCheck) ? "true \n\n" : "false \n\n");

        



    }
    for (int i = 0; i < solveForArray.size(); i++)
    {

        Js = solveForArray[i].J;
        Bxs = solveForArray[i].Bx;
        Bzs = solveForArray[i].Bz;
        performSolve(Js, Bxs, Bzs, &J, &Bx, &Bz, verbose);
    }

    MatDestroy(&J);
    MatDestroy(&Bx);
    MatDestroy(&Bz);
}