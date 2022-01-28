inline PetscInt getMatrixDim(PetscInt lattice_size)
{  
    return 2 << (lattice_size - 1);
}

inline PetscInt findLogOfTwo(PetscInt input)
{      
    PetscInt test = 1;
    counter = 1;  
    while (test < input) {test << 1; counter++;}
    return counter;
}


bool isPowerOfTwo(PetscInt input)
{
     return (!(input & (input - 1 )) && input);

}
void generateAdjacentBeta(PetscInt lattice_size, PetscInt i, PetscInt j, std::vector<double> &toAdd)
{
    PetscInt N = lattice_size * lattice_size;
    PetscInt matdim = getMatrixDim(lattice_size);
    std::vector<double> output;
    for (PetscInt i = 0; i < matdim; i++) 
    {
        output.push_back(0.0);
    }

    //q and q+1
    output[0] = 1.0;

    for (PetscInt x = 1; x < N + 1; x++){
        PetscInt twoRaisedX = 1;
        for (int i = 0; i < x; i++){
            twoRaisedX *= 2;
        }

        PetscInt thisStep = matdim / (twoRaisedX / 2);
        PetscInt nextStep = matdim / twoRaisedX;

        if (x == i || x == j){
            for (PetscInt i = 0; i < output.size(); i += thisStep){
                output[i + nextStep] = -1.0 * output[i];
            }
        }
        else{
            for (PetscInt i = 0; i < output.size(); i += thisStep){
                output[i + nextStep] = output[i];
            }
        }
    }
    bool writeMode = false;
    if (toAdd.size() < output.size())
        writeMode = true;
    for (PetscInt i = 0; i < output.size(); i++){
        if (writeMode)
            toAdd.push_back(output[i]);
        else
            toAdd.at(i) += output[i];
    }
}

void generateJMat(PetscInt lattice_size, Mat *matrix)
{
    PetscInt matrixDim = getMatrixDim(lattice_size);
    PetscInt N = lattice_size * lattice_size;
    std::vector<double> diagonal;
    std::vector<PetscInt> diagonal_map;
    for (PetscInt i = 0; i < matrixDim; i++){
        diagonal_map.push_back(i);
    }

    for (PetscInt q = 1; q < N + 1; q++){
        std::vector<double> addedDiagonal;

        PetscInt side_neighbor = q + 1;
        if(q % lattice_size == 0){
            side_neighbor -= lattice_size;
        }

        PetscInt down_neighbor = q + lattice_size;
        if(down_neighbor > N){
            down_neighbor = down_neighbor % N;
        }

        if(side_neighbor != q-1)
            generateAdjacentBeta(lattice_size, q, side_neighbor, addedDiagonal);

        if(down_neighbor != q - lattice_size)
            generateAdjacentBeta(lattice_size, q, down_neighbor, addedDiagonal);

        bool writeMode = addedDiagonal.size() > diagonal.size();
        for (PetscInt i = 0; i < addedDiagonal.size(); i++)
        {
            if (writeMode)
                diagonal.push_back(addedDiagonal[i]);
            else
                diagonal[i] += addedDiagonal[i];
        }
    }

    Vec v;

    VecCreate(PETSC_COMM_WORLD, &v);
    VecSetSizes(v, PETSC_DECIDE, matrixDim);
    VecSetFromOptions(v);
    VecSetValues(v, matrixDim, diagonal_map.data(), diagonal.data(), INSERT_VALUES);
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
    MatDiagonalSet(*matrix, v, INSERT_VALUES);
    VecDestroy(&v);
}

std::vector<double> generateBzPattern(PetscInt lattice_size, PetscInt level, double input)
{
    PetscInt N = lattice_size * lattice_size;
    if (level == N){
        std::vector<double> returner;
        returner.push_back(input);
        returner.push_back(input - 2);

        return returner;
    }
    else{
        std::vector<double> returner;

        std::vector<double> r1 = generateBzPattern(lattice_size, level + 1, input);
        std::vector<double> r2 = generateBzPattern(lattice_size, level + 1, input - 2);

        returner.insert(std::end(returner), std::begin(r1), std::end(r1));
        returner.insert(std::end(returner), std::begin(r2), std::end(r2));

        return returner;
    }
}

void generateBzMat(PetscInt lattice_size, Mat *matrix)
{
    PetscInt N = lattice_size * lattice_size;
    PetscInt matrixDim = getMatrixDim(lattice_size);

    //Levels of pattern recursion: lattice_size
    double start = N;
    std::vector<double> diagonal = generateBzPattern(lattice_size, 1, N);
    std::vector<PetscInt> diagonal_map;
    for (PetscInt i = 0; i < matrixDim; i++){
        diagonal_map.push_back(i);
    }
    Vec v;

    VecCreate(PETSC_COMM_WORLD, &v);
    VecSetSizes(v, PETSC_DECIDE, matrixDim);
    VecSetFromOptions(v);
    VecSetValues(v, matrixDim, diagonal_map.data(), diagonal.data(), INSERT_VALUES);
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
    MatDiagonalSet(*matrix, v, INSERT_VALUES);
    VecDestroy(&v);
}





void generateBxMat(PetscInt lattice_size, Mat *matrix)
{
    //Note
    //Row 0 of Bx: 1 only when column is a power of 2, including 1
    //Col 0 of Bx: 1 only when row is power of 2, including 1
    //0 aPetscInt diagonal
    //From Row 0: 1 aPetscInt diagonal where column is 1. There are as many ones as the column number. Then there are that many zeros
    PetscInt matrixDim = getMatrixDim(lattice_size);
    PetscInt N = lattice_size * lattice_size;


    PetscInt nonzeros = N * matrixDim;
    for (PetscInt i = 0; i < matrixDim; i++)
    {
        std::vector<PetscInt> colIndex;
        std::vector<double> row;
        
        for (PetscInt j = 0; j < i; j++){
            if (isPowerOfTwo(i - j )){
                //J determines if we insert a one or a zero.
                //One pattern: 1..0..1..
                PetscInt value = (((( j / (i-j)) + 1) % 2) == 0) ? 0 : 1;
                if(value == 1 ){
                    row.push_back(1);
                    colIndex.push_back(j);
                }
            }


        }
        
        for (PetscInt j = i; j < matrixDim; j++){
            if (isPowerOfTwo(j - i )){
                //J determines if we insert a one or a zero.
                //One pattern: 1..0..1..
                PetscInt value = (((( i / (j-i)) + 1) % 2) == 0) ? 0 : 1;
               
                if(value == 1 ){
                    row.push_back(1);
                    colIndex.push_back(j);
                }
            }
        }
        PetscInt rows_to_set[] = {i};
        MatSetValues(*matrix, 1, &rows_to_set[0], colIndex.size(), colIndex.data(), row.data(), INSERT_VALUES);
            
    }
}
