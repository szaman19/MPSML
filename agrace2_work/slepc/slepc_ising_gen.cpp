#include<stdio.h>

#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"
#include<iostream>

static char help[] = "Help message";

//Constant Pauli matrices
Mat SigmaXI;
Mat SigmaZI;
Mat ID;

//Create the Pauli Matrices
void createSigmaXI(Mat * SigmaXI){
    MatCreate(PETSC_COMM_WORLD, SigmaXI);
    MatSetSizes(*SigmaXI,PETSC_DECIDE, PETSC_DECIDE, 2, 2 );
    MatSetFromOptions(*SigmaXI);
    //Set Values
    MatSetValue(*SigmaXI, 0,0, 0.0, INSERT_VALUES);
    MatSetValue(*SigmaXI, 1,0, 1.0, INSERT_VALUES);
    MatSetValue(*SigmaXI, 0,1, 1.0, INSERT_VALUES);
    MatSetValue(*SigmaXI, 1,1, 0.0, INSERT_VALUES);
    MatAssemblyBegin(*SigmaXI);
    MatAssemblyEnd(*SigmaXI);
}
void createID(Mat * ID){
    MatCreate(PETSC_COMM_WORLD, ID);
    MatSetSizes(*ID,PETSC_DECIDE, PETSC_DECIDE, 2, 2 );
    MatSetFromOptions(*ID);
    //Set Values
    MatSetValue(*ID, 0,0, 1.0, INSERT_VALUES);
    MatSetValue(*ID, 1,0, 0.0, INSERT_VALUES);
    MatSetValue(*ID, 0,1, 0.0, INSERT_VALUES);
    MatSetValue(*ID, 1,1, 1.0, INSERT_VALUES);
    MatAssemblyBegin(*ID);
    MatAssemblyEnd(*ID);
}

void createSigmaZI(Mat * SigmaZI){
    MatCreate(PETSC_COMM_WORLD, SigmaZI);
    MatSetSizes(*SigmaXZ,PETSC_DECIDE, PETSC_DECIDE, 2, 2 );
    MatSetFromOptions(*SigmaXI);
    //Set Values
    MatSetValue(*SigmaZI, 0,0, 1.0, INSERT_VALUES);
    MatSetValue(*SigmaZI, 1,0, 0.0, INSERT_VALUES);
    MatSetValue(*SigmaZI, 0,1, 0.0, INSERT_VALUES);
    MatSetValue(*SigmaZI, 1,1, -1.0, INSERT_VALUES);
    MatAssemblyBegin(*SigmaZI);
    MatAssemblyEnd(*SigmaZI);
}

//Tensor
//Finds the Kronecker Tensor Product with 3 matrices
//A- first one
//B- second one
//C- output
//Requires that A and b be assembled, C be unassembled and not created
void tensor(Mat * A, Mat * B, Mat * C, int aRows, int aCols, int bRows, int bCols){
    MatCreate(PETSC_COMM_WORLD, C);
    MatSetSizes(*C, PETSC_DECIDE, PETSC_DECIDE, aRows * bRows, aCols * bCols);
    MatSetFromOptions(*C);
    for(int i = 0; i < aRows; i++){
        for(int j = 0; j < aCols; j++){
            double multiplier[1];
            MatGetValues(A, 1, 1, i, j, &multiplier);
            for(int s = 0; s < bRows; s++){
                for(int t = 0; t < bCols; t++){
                    double bValue[1];
                    MatGetValues(B, 1, 1, s, t, &bValue);
                    MatSetValue(*C, i * bRows + s, j * bCols + t, multiplier[0] * bValue[0] );
                }
            }
        }
    }
}

//Helper integer function
int intPower(int input , int power){
    int out = 1;
    for(int x = 0; x < power; x++){
        out *= input;
    }
    return out;
}

//Generate self interaction
//Helper function for generating Bx and Bz matrices
//Changes the matrix output, which should be created but not set.
//Sets its size but does NOT assemble it
void generateSelfInteraction(Mat* output, char spin, int i, int N){
    //calculate size of final matrix
    int finalX = intPower(2,N);
    int finalY = intPower(2,N);

    //use these as storage so we don't need to muck about with matrices every iteration
    int currentX = 1;
    int currentY = 1;

    //Set up output
    MatSetSizes(*output, PETSC_DECIDE, PETSC_DECIDE, finalX, finalY);
    MatSetFromOptions(*output);
    MatSetValue(*output, 0, 0, 1.0);

    for(int x = 1; x < N + 1; x++){
        //Now we need to juggle Ouput and assemble it so it works with the tensor function as A
        Mat copyOfOutput;
        MatCopy(output, copyOfOutput, SAME_NONZERO_PATTERN);
        MatAssemblyBegin(copyOfOutput);
        MatAssemblyEnd(copyOfOutput);

        //Create the temp output but don't assemble it
        Mat tensorProductOutput;
        if(x == i){
            if(spin == 'x'){
                tensor(output, SigmaXI, tensorProductOutput, currentSize, currentSize, 2, 2);
            }
            else{
                tensor(output, SigmaZI, tensorProductOutput, currentSize, currentSize, 2, 2);
            }
        }
        else{
            tensor(output, ID, tensorProductOutput, currentSizeX, currentY,2,2);
        }

        currentSizeX *=2;
        currentSizeY *=2;
        MatAssemblyBegin(tensorProductOutput);
        MatAssemblyEnd(tensorProductOutput);
        //Zero out the output matrix, copy whats in copy to the output, destroy copy
        MatZeroEntries(output);
        
        for(int x = 0; x< currentSizeX; x++){
            for(int y = 0; y < currentSizeY; y++){
                double bValue[1];
                MatGetValues(tensorProductOutput, 1, 1, x, y, &bValue);
                MatSetValue(output, x, y, bValue[0], INSERT_VALUES);
            }
        }
        MatDestroy(tensorProductOutput);
        MatDestroy(copyOfOutput);
    }

}

//Generate self interaction
//Helper function for generating J matrices
//Changes the matrix output, which should be created but not set.
//Sets its size but does NOT assemble it
void generateAdjacentInteraction(Mat* output, int i, int j, int N){
    //calculate size of final matrix
    int finalX = intPower(2,N);
    int finalY = intPower(2,N);

    //use these as storage so we don't need to muck about with matrices every iteration
    int currentX = 1;
    int currentY = 1;

    //Set up output
    MatSetSizes(*output, PETSC_DECIDE, PETSC_DECIDE, finalX, finalY);
    MatSetFromOptions(*output);
    MatSetValue(*output, 0, 0, 1.0);

    for(int x = 1; x < N + 1; x++){
        //Now we need to juggle Ouput and assemble it so it works with the tensor function as A
        Mat copyOfOutput;
        MatCopy(output, copyOfOutput, SAME_NONZERO_PATTERN);
        MatAssemblyBegin(copyOfOutput);
        MatAssemblyEnd(copyOfOutput);

        //Create the temp output, but don't assemble it
        Mat tensorProductOutput;
        if(x == i || x == j){
            tensor(output, SigmaZI, tensorProductOutput, currentSize, currentSize, 2, 2);
        }
        else{
            tensor(output, ID, tensorProductOutput, currentSizeX, currentY,2,2);
        }

        currentSizeX *=2;
        currentSizeY *=2;
        MatAssemblyBegin(tensorProductOutput);
        MatAssemblyEnd(tensorProductOutput);
        //Zero out the output matrix, copy whats in copy to the output, destroy copy
        MatZeroEntries(output);
        
        for(int x = 0; x< currentSizeX; x++){
            for(int y = 0; y < currentSizeY; y++){
                double bValue[1];
                MatGetValues(tensorProductOutput, 1, 1, x, y, &bValue);
                MatSetValue(output, x, y, bValue[0], INSERT_VALUES);
            }
        }
        MatDestroy(tensorProductOutput);
        MatDestroy(copyOfOutput);
    }


}

//Generates the J Matrix
//Expects JTerms to be created, but not set up
void generateJMatrix(Mat * JTermsOutput, int N, int latice_size_one_dimension){
    int matrixDim = intPower(2 * N)
    MatSetSizes(*JTermsOutput, PETSC_DECIDE, PETSC_DECIDE, matrixDim, matrixDim);
    MatSetFromOptions(*JTermsOutput);
    MatAssemblyBegin(*JTermsOutput);
    MatAssemblyEnd(*JTermsOutput);
    
    for (int q = 1; q < N+1; q++)
    {
        if (q + 1 < (N+1) && (q + 2) % latice_size_one_dimension != 0){
            //JTerms = JTerms + generateAdjacentInteractionZ(q, q+1);
            Mat tempAdjInteraction;
            MatCreate(PETSC_COMM_WORLD, &tempAdjInteraction);
            generateAdjacentInteraction(&tempAdjInteraction, q, q+1, N);
            MatAssemblyBegin(tempAdjInteraction);
            MatAssemblyEnd(tempAdjInteraction);
            MatAXPY(*JTermsOutput, 1, tempAdjINteraction, SAME_NONZER_PATTERN);
            MatDestroy(tempAdjInteraction);
        }
        if (q + latice_size_one_dimension < (N+1)){
            // JTerms = JTerms + generateAdjacentInteractionZ(q, q + latice_size_one_dimension);
            Mat tempAdjInteraction;
            MatCreate(PETSC_COMM_WORLD, &tempAdjInteraction);
            generateAdjacentInteraction(&tempAdjInteraction, q, q+latice_size_one_dimesion, N);
            MatAssemblyBegin(tempAdjInteraction);
            MatAssemblyEnd(tempAdjInteraction);
            MatAXPY(*JTermsOutput, 1, tempAdjINteraction, SAME_NONZER_PATTERN);
            MatDestroy(tempAdjInteraction);
        }
           
    }
}

//Generates the Bx Matrix
//Expects BxTermsOut to be created, but not set up
void generateBxMatrix(Mat * BxTermsOut, int N, int latice_size_one_dimension){
    int matrixDim = intPower(2 * N)
    MatSetSizes(*BxTermsOut, PETSC_DECIDE, PETSC_DECIDE, matrixDim, matrixDim);
    MatSetFromOptions(*BxTermsOut);
    MatAssemblyBegin(*BxTermsOut);
    MatAssemblyEnd(*BxTermsOut);
    for (int q = 1; q < N+1; q++)
    {
        Mat tempSelfInteraction;
        MatCreate(PETSC_COMM_WORLD, &tempSelfInteraction);
        generateSelfInteraction(&tempSelfInteraction, 'x', q, int N);
        MatAssemblyBegin(tempSelfInteraction);
        MatAssemblyEnd(tempSelfInteraction);
        MatAXPY(*JTermsOutput, 1, tempSelfInteraction, SAME_NONZER_PATTERN);
        MatDestroy(tempSelfInteraction);
    }


}

//Generates the Bz Matrix
//Expects BzTermsOut to be created, but not set up
void generateBzMatrix(Mat * BzTermsOut, int N, int latice_size_one_dimension){
    int matrixDim = intPower(2 * N)
    MatSetSizes(*BzTermsOut, PETSC_DECIDE, PETSC_DECIDE, matrixDim, matrixDim);
    MatSetFromOptions(*BzTermsOut);
    
    for (int q = 1; q < N+1; q++)
    {
        Mat tempSelfInteraction;
        MatCreate(PETSC_COMM_WORLD, &tempSelfInteraction);
        generateSelfInteraction(&tempSelfInteraction, 'z', q, int N);
        MatAssemblyBegin(tempSelfInteraction);
        MatAssemblyEnd(tempSelfInteraction);
        MatAXPY(*JTermsOutput, 1, tempSelfInteraction, SAME_NONZER_PATTERN);
        MatDestroy(tempSelfInteraction);
    }

    
}

int main(int argc, char * argv[]){

    if(argc != 5){
        std::cout << "Usage: ./ <executable name goes here> <latice side size> < J value > < Bsubx Value> < Bsubz value > "<< std::endl;
        return 0;
    }

    // Hamiltonian generation Parameters
    int latice_size = 2;
    double J = 1.0;
    double Bx = 1;
    double Bz = 1;
    /*
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
    */

    std::cout << "Initializing Slepc" <<std::endl;
    SlepcInitialize(&argc, &argv, (char*) 0, help);

    //Initialize the Pauli Matrices
    createSigmaXI(&SigmaXI);
    createSigmaZI(&SigmaZI);
    createID(&ID);
    std::cout << "Initialized Pauli Matrices" <<std::endl;
    
    Mat BzTerms;
    Mat BxTerms;
    Mat JTerms;
    MatCreate(PETSC_COMM_WORLD, BzTerms);
    MatCreate(PETSC_COMM_WORLD, BxTerms);
    MatCreate(PETSC_COMM_WORLD, JTerms);

    generateJMatrix(&JTerms, , latice_size)


    EPS eps;                // Eigensolver context
    Mat A;                  // The matrix that will eventually hold the hamiltonian. 
    Vec xr, xi;             // Eigenvectors x
    PetscScalar kr, ki, im, re;     // Eigenvalues, k
    PetscInt j, nconv;      // Petsc parameters
    PetscReal error;        // Errorvalue

    
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
