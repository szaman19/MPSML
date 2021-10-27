#include "SlepcMat.h"


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
        MatGetValues(mat, 1, r, 1, c, &existingValue[0]);
        //New Value in Matrix = Existing Value + (Value - Existing Value)
        value -= existingValue[0];
        Mat temp;
        MatCreate(PETSC_COMM_WORLD, &temp);
        MatSetSizes(temp, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
        MatSetFromOptions(temp);
        MatSetUp(temp);
        int ownershiprow1;
        int ownershiprow2;
        MatGetOwnershipRange(temp,&ownershiprow1,&ownershiprow2);
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

    for(PetscInt i = 0; i < oldMatrix.rows; i++){
        for(PetscInt j = 0; j < oldMatrix.cols; j++){
            //Get the multiplier
            double multiplier[1];
            PetscInt rowToGetM[] = {i};
            PetscInt colToGetM[] = {j};
            
        
            MatGetValues(oldMatrix.mat, 1, rowToGetM, 1, colToGetM, &multiplier[0]);
                for(PetscInt x = 0; x < other.rows; x++){
                    for(PetscInt y = 0; y < other.cols; y++){
                        double otherMatrixValue[1];
                        PetscInt rowToGet[] = {x};
                        PetscInt colToGet[] = {y};

                        int ownOther1;
                        int ownOther2;
                        MatGetOwnershipRange(other.mat,&ownOther1,&ownOther2);
                        
                        MatGetValues(other.mat, 1, rowToGet, 1, colToGet, &otherMatrixValue[0]);
                        //We can't use the built in method I wrote because mat is not yet assembled, so we use this instead
                        MatSetValue(mat, i * other.rows + x, j * other.cols + y, multiplier[0] * otherMatrixValue[0], INSERT_VALUES);
                        
                    }
                }
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
    /* Put the results into output */
    std::cout << nconv << std::endl;
    EPSGetEigenpair(solver, 0, vaReal, vaImaginary, *veReal, *veImaginary);
    EPSDestroy(&solver);
}