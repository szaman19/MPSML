#include<algorithm>

bool checkHermitian(Mat * test){
    PetscBool hermit;
    MatIsHermitian(*test, 0.0, &hermit);
    return hermit == PETSC_TRUE;
}

PetscInt findPowerOfTwo(PetscInt raisedTo){
    if(raisedTo < 0 ) raisedTo *=-1;
    PetscInt rv = 1;
    for(PetscInt i = 0; i < raisedTo; i++){
        rv *=2;
    }
    return rv;
}

bool checkAgainstShehtabsMath(Mat * test, PetscInt lattice_size){
    PetscInt matrixDim = getMatrixDim(lattice_size);
    PetscInt N = lattice_size * lattice_size;

    Mat Z;
    MatCreate(PETSC_COMM_WORLD, &Z);
    MatSetSizes(Z, PETSC_DECIDE, PETSC_DECIDE, getMatrixDim(lattice_size), getMatrixDim(lattice_size));
    MatSetFromOptions(Z);
    MatSetUp(Z);

    for(PetscInt i = 0; i < matrixDim; i++){
        std::vector<PetscInt> cols;
        std::vector<double> values;
        for(PetscInt h = 0; h < N; h++){
            for(PetscInt j = 0; j < matrixDim; j++){
                if( ((i == findPowerOfTwo( N - h - 1) + j ) && j < findPowerOfTwo(N - h - 1)) 
                ||  ((j == findPowerOfTwo( N - h - 1) + i ) && i < findPowerOfTwo(N - h - 1)) ){
                    if(std::find(cols.begin(), cols.end(), j) == std::end(cols)){
                        cols.push_back(j);
                        values.push_back(1.0);
                    }
                }
            }

        }
        PetscInt rows_to_set[] = {i};
        MatSetValues(Z, 1, &rows_to_set[0], cols.size(), cols.data(), values.data(), INSERT_VALUES);

    }

    
    MatAssemblyBegin(Z, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Z, MAT_FINAL_ASSEMBLY);
    MatView(Z, PETSC_VIEWER_STDOUT_WORLD );
    MatView(*test, PETSC_VIEWER_STDOUT_WORLD );
    PetscBool result;
    MatIsTranspose(*test, Z, 0.0, &result);
    MatDestroy(&Z);
    return result == PETSC_TRUE;




}