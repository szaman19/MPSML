/* Objective-oriented SLEPc and PETSc Wrapper for C++ */
/* Andrew Grace, 2021                                 */

#include<iostream>
#include<stdio.h>
#include"slepcsys.h"
#include<slepceps.h>
#include"petscsys.h"


class SlepcMat{
    private:
        Mat mat;

        /*  
        PETSC matrices have several stages to getting an actual usable matrix.
        Step 0: Declaring the variable. Should be done automatically by this class
        Step 1: Creating the matrix with MatCreate. This allocates space
        Step 2: Setting the sizes
        Step 3: Setting the matrix. Corresponds to MatSetFromOptions
        Step 4: Setting the values with MatSetValue
        Step 5: Matrix Assembly

        My wrapper uses 4 booleans to keep track of each step. 
        */
        bool isCreated;
        bool isSizeSet;
        bool isOptionsSet;
        bool isAssembled;

        /*
        It also stores some other data about the matrix, including row and col count
        */

        int rows;
        int cols;
        
        /*
        Basic translators for PETSC
        */
        
        
        
        void assembleMatrix(){
            MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
            isAssembled = true;
        }

        /*
        Disassemble is useful for other operations besides add and subtract
        */
        void destroy(){
            MatDestroy(&mat);
            isCreated = false;
            isSizeSet = false;
            isOptionsSet = false;
            isAssembled = false;
            rows = 0;
            cols = 0;
        }

        /*
        Autocomplete function
        Performs assembly if required,
        Sets matrix to zero matrix if isOptionsSet is not true
        */

        void autocomplete(){
            bool setToZero = false;
            if(!isCreated){
                isSizeSet = false;
                isOptionsSet = false;
                isAssembled = false;
                createMatrix();
            }
            if(!isSizeSet){
                setSizes(rows, cols);
                isOptionsSet = false;
                isAssembled = false;
            }
            if(!isOptionsSet){
                setOptions();
                setToZero = true;
            }
            if(setToZero) MatZeroEntries(mat);
            if(!isAssembled){
                assembleMatrix();
            }
        }

        int getMatBegin(){
            int begin;
            int end;
            MatGetOwnershipRange(mat, &begin, &end);
            return begin;
        }

        int getMatEnd(){
            int begin;
            int end;
            MatGetOwnershipRange(mat, &begin, &end);
            return begin;
        }

    public:
        SlepcMat(){
            isCreated = false;
            isSizeSet = false;
            isOptionsSet = false;
            isAssembled = false;
            rows = 0;
            cols = 0;
        }

        SlepcMat(int rows, int cols){
            isCreated = false;
            isSizeSet = false;
            isOptionsSet = false;
            isAssembled = false;
            createMatrix();
            setSizes(rows, cols);
            //Can't do anything else
        }

        SlepcMat(const SlepcMat &other){
            this->rows = other.rows;
            this->cols = other.cols;
            createMatrix();
            setSizes(this->rows, this->cols);
            setOptions();
            assembleMatrix();
            MatConvert(other.mat, MATSAME, MAT_INITIAL_MATRIX, &this->mat);
            assembleMatrix();
        }

        ~SlepcMat(){
            destroy();
        }

        void setValue(double value, int row, int col){
            if(!isOptionsSet){
                setOptions();
            }

            if(!isAssembled && isSizeSet && isCreated && isOptionsSet){
                int begin = getMatBegin();
                MatSetValue(mat, begin + row, col, value, INSERT_VALUES);
               
            }
            else {
                /*
                SlepcMat burner(*this);
                destroy();
                MatSetValue(burner.mat, row, col, value, INSERT_VALUES);
                createMatrix();
                setSizes(this->rows, this->cols);
                setOptions();
                MatDuplicate(burner.mat, MAT_COPY_VALUES, &mat);
                */
            }
            
        }

        void setTo ( SlepcMat other){
            destroy();
            this->rows = other.rows;
            this->cols = other.cols;
            other.autocomplete();
            createMatrix();
            setSizes(this->rows, this->cols);
            setOptions();
            
            MatConvert(other.mat, MATSAME, MAT_INITIAL_MATRIX, &this->mat);
            assembleMatrix();
        }

        SlepcMat& operator= (const SlepcMat& other){
            if(this != &other){
                this->rows = other.rows;
                this->cols = other.cols;
                createMatrix();
                setSizes(this->rows, this->cols);
                setOptions();
                MatConvert(other.mat, MATSAME, MAT_INITIAL_MATRIX, &this->mat);
            }
            return *this;
        }

        SlepcMat operator+ (SlepcMat other){
            
            SlepcMat output(*this);
            if(this->rows = other.rows && this->cols == other.cols){
                other.autocomplete();
                autocomplete();
                MatAXPY(output.mat, 1.0, other.mat, SAME_NONZERO_PATTERN);
            }
            return output;
        }

        void createMatrix(){
            MatCreate(PETSC_COMM_WORLD, &mat);
            isCreated = true;
        }

        void setSizes(int rows, int cols){
            MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
            this->rows = rows;
            this->cols = cols;
            isSizeSet = true;
        }
        void setOptions(){
            MatSetFromOptions(mat);
            MatSetUp(mat);
            isOptionsSet = true;
        }

        

        SlepcMat operator* (double x){
            SlepcMat output(*this);
            SlepcMat zeroMat(this->rows, this->cols);
            zeroMat.createMatrix();
            zeroMat.setSizes(this->rows, this->cols);
            zeroMat.setOptions();
            MatZeroEntries(zeroMat.mat);
            this->autocomplete();
            zeroMat.autocomplete();
            MatAYPX(output.mat, x, zeroMat.mat, SAME_NONZERO_PATTERN);
            return output;
        }

        void tensor_destroy_previous_matrix( SlepcMat & other){
            //Make a copy of this matrix.
            autocomplete();
            SlepcMat oldMatrix(*this);

            other.autocomplete();
            //Destroy the old matrix
            destroy();
            //Compute size of new matrix
            rows = oldMatrix.rows * other.rows;
            cols = oldMatrix.cols * other.cols;
            createMatrix();
            setSizes(rows, cols);
            setOptions();
            //Compute the tensor product
            for(PetscInt i = 0; i < oldMatrix.rows; i++){
                for(PetscInt j = 0; j < oldMatrix.cols; j++){
                    double multiplier[1];
                    PetscInt rowToGetM[] = {i};
                    PetscInt colToGetM[] = {j};
                    MatGetValues(oldMatrix.mat, 1, rowToGetM, 1, colToGetM, &multiplier[0]);

                    for(PetscInt x = 0; x < other.rows; x++){
                        for(PetscInt y = 0; y < other.cols; y++){
                            double otherMatrixValue[1];
                            PetscInt rowToGet[] = {x};
                            PetscInt colToGet[] = {y};
                            MatGetValues(other.mat, 1, rowToGet, 1, colToGet, &otherMatrixValue[0]);
                            MatSetValue(mat, i * other.rows + x, j * other.cols + y, multiplier[0] * otherMatrixValue[0], INSERT_VALUES);
                            }
                    }
                }
            }
            assembleMatrix();
            oldMatrix.destroy();
        }
        

        void view(){
            if(!isAssembled) this->autocomplete();
            MatView(mat, PETSC_VIEWER_STDOUT_WORLD);
        }

        void getLowestEigenpair(PetscScalar * vaReal, PetscScalar * vaImaginary, Vec * veReal, Vec * veImaginary){
            /* Assemble the matrix if we haven't already */
            autocomplete();
            
            /* Set up SLEPC */
            EPS solver;
            int nconv;
            EPSCreate(PETSC_COMM_WORLD, &solver);
            EPSSetOperators(solver, mat, NULL);
            EPSSetProblemType(solver, EPS_NHEP);
            EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL);
            EPSSetFromOptions(solver);
            EPSSolve(solver);
            EPSGetConverged(solver, &nconv);

            /* Put the results into output */

            EPSGetEigenpair(solver, 0, vaReal, vaImaginary, veReal, veImaginary);

            EPSDestroy(&solver);


        }


        





};

