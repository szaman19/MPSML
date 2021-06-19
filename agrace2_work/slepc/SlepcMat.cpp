#include "SlepcMat.h"
#include<cmath>
#include<limits>
#include<mpi.h>
#include <unistd.h>

/* Create the matrix by running the standard PETSC routine for making a matrix. Store it in mat */
void SlepcMat::createMatrix(){
    MatCreate(PETSC_COMM_WORLD, &mat);
    MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
    MatSetFromOptions(mat);
    MatSetUp(mat);
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    assembledState = true;
}

/* Destroy the matrix */
void SlepcMat::destroy(){
    assembledState = false;
    MatDestroy(&mat);
    rows = 1;
    cols = 1;
}

/* No-value constructor makes a 1x1 matrix */
SlepcMat::SlepcMat(){
    rows = 1;
    cols = 1;
    createMatrix();
}

/* RowxCols constructor */
SlepcMat::SlepcMat(int rows, int cols){
    this->rows = rows;
    this->cols = cols;
    createMatrix();

}

/* Copy constructor */
SlepcMat::SlepcMat(const SlepcMat &other){

    this->rows = other.rows;
    this->cols = other.cols;
    createMatrix();
    MatDuplicate(other.mat, MAT_COPY_VALUES, &this->mat);
}

/* Destructor */
SlepcMat::~SlepcMat(){
    
    if(assembledState){
        destroy();

    }
    else{

    }
    
}

/* Set a value in the matrix. First get the value in the spot, then calculate the difference between the chosen new value and the old value. */
/* Using math, we get newValue = (newValue-oldValue) + oldValue */
/* Make a separate matrix with that value in the correct spot and add it so we don't have to dissassemble the matrix */
void SlepcMat::setValue(double value, int row, int col){
    if(row < rows && col < cols){
        double existingValue[1];
        PetscInt c[] = {col};
        PetscInt r[] = {row};
        int ownershiprow1;
        int ownershiprow2;
        MatGetOwnershipRange(mat,&ownershiprow1,&ownershiprow2);
        if(row >= ownershiprow1 && row < ownershiprow2){
            MatGetValues(mat, 1, r, 1, c, &existingValue[0]);
        }
        //New Value in Matrix = Existing Value + (Value - Existing Value)
        value -= existingValue[0];
        Mat temp;
        MatCreate(PETSC_COMM_WORLD, &temp);
        MatSetSizes(temp, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
        MatSetFromOptions(temp);
        MatSetUp(temp);
        
        if(row >= ownershiprow1 && row < ownershiprow2){
            MatSetValue(temp, row, col, value, INSERT_VALUES);
        }
        MatAssemblyBegin(temp, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(temp, MAT_FINAL_ASSEMBLY);
        MatAYPX(mat, 1.0, temp, DIFFERENT_NONZERO_PATTERN);
        MatDestroy(&temp);
    }
}

SlepcMat& SlepcMat::operator= (const SlepcMat& other){
    if(this != &other){
        this->rows = other.rows;
        this->cols = other.cols;
        createMatrix();
        MatDuplicate(other.mat, MAT_COPY_VALUES, &this->mat);
    }
    return *this;
}

SlepcMat SlepcMat::operator+ (SlepcMat other){

    SlepcMat output(*this);

    if(this->rows = other.rows && this->cols == other.cols){

        MatAXPY(output.mat, 1.0, other.mat, DIFFERENT_NONZERO_PATTERN);

    }
    return output;
}

SlepcMat SlepcMat::operator* (double x){
    SlepcMat output(*this);
    SlepcMat zeroMat(this->rows, this->cols);
    zeroMat.createMatrix();
    MatZeroEntries(zeroMat.mat);
    MatAYPX(output.mat, x, zeroMat.mat, SAME_NONZERO_PATTERN);
    return output;
}

void SlepcMat::tensor_destroy_previous_matrix( SlepcMat & other){
    //Make a copy of this matrix.
    SlepcMat oldMatrix(*this);
    //Destroy the old matrix
    destroy();
    //Compute size of new matrix and set it up
    rows = oldMatrix.rows * other.rows;
    cols = oldMatrix.cols * other.cols;
    MatCreate(PETSC_COMM_WORLD, &mat);
    MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
    MatSetFromOptions(mat);
    MatSetUp(mat);
    //Compute the tensor product
    //PETSC does not like it when you try to do this routine over multiple nodes!

    //Let's try adding some MPI code to fix that.
    //Local copy of the tensor factor.
    double tensorFactor[other.rows][other.cols];
    //Populate the copy
    //std::cout <<  "made it to shareMatrix" << std::endl;
    shareMatrix(&tensorFactor[0][0], &other.mat, other.rows, other.cols);
    //std::cout <<  "made it past shareMatrix" << std::endl;
    //Get the ownership range of the old matrix.
    int ownOldBegin;
    int ownOldEnd;
    MatGetOwnershipRange(oldMatrix.mat,&ownOldBegin,&ownOldEnd);
    
    for(PetscInt i = 0; i < oldMatrix.rows; i++){
        for(PetscInt j = 0; j < oldMatrix.cols; j++){
            //Do we own this value???
            if(i >= ownOldBegin && i < ownOldEnd){
                //If we do, get this value as the multiplier
                double multiplier[1];
                PetscInt rowToGet[] = {i};                        
                PetscInt colToGet[] = {j};
                //std::cout << "attempting to get multiplier where i=" << i << ", j=" << j << std::endl;
                MatGetValues(oldMatrix.mat, 1, rowToGet, 1, colToGet, &multiplier[0]);
                //PetscPrintf(MPI_COMM_WORLD, "array checks: i=%d, j=%d\n", rowToGet[0], colToGet[0] );   
                //And then use our local copy to complete the tensor product
                
                /*
                for(PetscInt x = 0; x < other.rows; x++){
                    //PetscSynchronizedFPrintf(MPI_COMM_WORLD, stderr, "attempting to set things\n"); 
                    for(PetscInt y = 0; y < other.cols; y++){
                       // PetscSynchronizedFPrintf(MPI_COMM_WORLD, stderr, "attempting to set things\n"); 
                        //sleep(1);
                        //PetscSynchronizedFPrintf(MPI_COMM_WORLD, stderr, "attempting to set x=%d, y=%d\n", x,  y);   
                        
                        MatSetValue(mat, i * other.rows + x, j * other.cols + y, multiplier[0] * tensorFactor[x][y], INSERT_VALUES);
                    }
                }
                */
                //Rebuilt to use MatSetValues
                //Assemble the array of rows to be updated and multiply tensorFactor by the multiplier
                
                PetscInt xVals[other.rows];
                PetscInt yVals[other.cols];
                double * tempMultiplied = new double[other.rows * other.cols];
                for(int x = 0; x < other.rows; x++){
                    xVals[x] = i * other.rows + x;
                    for(int y = 0; y < other.cols; y++){
                        if(x == 0){
                            yVals[y] = j * other.cols + y;
                        }
                        tempMultiplied[other.cols * x + y] = multiplier[0] * tensorFactor[x][y];
                    }
                }
                //MatSetValues might be faster.
                MatSetValues(mat, other.rows, xVals, other.cols, yVals, (const double *) tempMultiplied, INSERT_VALUES);
                
            }
            
            //If we don't, then don't do anything.      
        }
    }
    
    //Clean up and assemble
    oldMatrix.destroy();
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
}

void SlepcMat::view(){
    MatView(mat, PETSC_VIEWER_STDOUT_WORLD);
}

void SlepcMat::getLowestEigenpair(PetscScalar * vaReal, PetscScalar * vaImaginary, Vec * veReal, Vec * veImaginary){
    
    /* Set up SLEPC */
    EPS solver;
    int nconv;
    /* Initialize vectors */
    MatCreateVecs(mat, NULL, veImaginary);
    MatCreateVecs(mat, NULL, veReal);

    EPSCreate(PETSC_COMM_WORLD, &solver);
    EPSSetOperators(solver, mat, NULL);
    EPSSetProblemType(solver, EPS_NHEP);
    EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL);
    EPSSetFromOptions(solver);
    EPSSolve(solver);
    EPSGetConverged(solver, &nconv);
    EPSGetEigenpair(solver, 0, vaReal, vaImaginary, *veReal, *veImaginary);
    EPSDestroy(&solver);
}



void SlepcMat::shareMatrix(double * matrixArray, Mat * matrix, int numRows, int numCols){
    //Use MPI to get values in a SLEPC matrix and created copies on each node in matrixArray
    
    int numberOfNodes;
    int thisNode;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisNode);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);

    //See what parts of the matrix this processor owns
    int matrixRowOwnershipBegin;
    int matrixRowOwnershipEnd;
    MatGetOwnershipRange(*matrix,&matrixRowOwnershipBegin,&matrixRowOwnershipEnd);

    //Use MPI_AllGather to populate the shared array.
    //Credit to Shehtab Zaman, for finding that solution.

    //Temporary matrix buffer required.
    //PETSC should give each process an equal number of rows
    //Using that fact, we have that the temporary buffer should have approximately numRows / numberOfNodes rows and numCols cols

    int tempRowSize = numRows / numberOfNodes;
    if(numRows % numberOfNodes != 0) tempRowSize++;
    if(tempRowSize == 0){
        tempRowSize++;
        //printf("tempRowSize = 1\n");
    }
    double * tempBuffer = new double[tempRowSize * numCols];
    

    //Also create a buffer that can hold numberOfNodes * tempRowSize nodes
    double* tempAllBuffer = new double[(tempRowSize * numberOfNodes) * numCols];

    //Use numRowsProcessed as the row index for the matrix
    int numRowsProcessed = 0;
    for(int i = 0; i < numRows && numRowsProcessed < tempRowSize; i++){
        if(i < matrixRowOwnershipBegin){
            //Do nothing if we havent reached the range yet
        }
        else if(i >= matrixRowOwnershipBegin && i < matrixRowOwnershipEnd){
            //Get a value from the matrix and place it into our first (local) temporary buffer
            for(int j = 0; j < numCols; j++){
                PetscInt rowToGet[] = {i};
                PetscInt colToGet[] = {j};
                double value[1];
                MatGetValues(*matrix, 1, rowToGet, 1, colToGet, &value[0]);
                tempBuffer[numCols * numRowsProcessed + j] = value[0];
            }
            numRowsProcessed++;
        }

        else if(i >= matrixRowOwnershipEnd){
            //If we ever reach this, then temprowsize is 1 + the actual number of rows owned
            //We signal this by setting the first value to NAN.
            tempBuffer[numRowsProcessed * numCols] = std::numeric_limits<double>::quiet_NaN();
            numRowsProcessed++;
        }
    }

    //The previous approach went through all the nodes and every matrix value, creating local copies with MPI_Bcast
    //AllGather reduces the complexity of this operation.
    //Sync up, copy, sync again, and ditch the temp buffer

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(tempBuffer, tempRowSize * numCols, MPI_DOUBLE, tempAllBuffer, tempRowSize * numCols, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] tempBuffer;

    //Clean the data of NAN rows; if a row's first value is NaN, then this row was a gap row. 
    int currentNonNANRowCount = 0;

    for(int i = 0; i < tempRowSize; i++){
        bool isAllNan = std::isnan(tempAllBuffer[numCols * i]);
          
        for(int j = 0; j < numCols && isAllNan == false; j++){
            //std::cout << "node " << thisNode << " setting (" << currentNonNANRowCount << ", " << j << ") as " << tempAllBuffer[numCols * i + j] << std::endl;
            matrixArray[numCols * currentNonNANRowCount + j] = tempAllBuffer[numCols * i + j];
        }
        currentNonNANRowCount += (isAllNan) ? 0 : 1;
    }

    delete[] tempAllBuffer;
    MPI_Barrier(MPI_COMM_WORLD);

}