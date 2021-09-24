#include<stdio.h>
#include<Eigen/sparse>
#include<vector>
#include"DynamicMatrix.h"

typedef Eigen::SparseMatrix<double> SpMat;

std::vector<double> generateBzPattern(int level, long N, double input){

    if(level == N){
        std::vector <double> returner;
        returner.push_back(input);
        returner.push_back(input - 2);
        return returner;
    }
    else{
        std::vector<double> returner;
        std::vector<double> r1 = generateBzPattern(level + 1, N, input);
        std::vector<double> r2 = generateBzPattern(level + 1, N, input - 2);
        returner.insert(std::end(returner), std::begin(r1), std::end(r1));
        returner.insert(std::end(returner), std::begin(r2), std::end(r2));
        return returner;
    }
}

void generateBzMatrix(long lattice_size, Eigen::SparseMatrix<double> *mat ){
    //Get the values and put them into triplets
    long N = lattice_size * lattice_size;
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    //Levels of pattern recursion: lattice_size
    double start = N;
    std::vector<double> diagonals = generateBzPattern(1,N, N);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < diagonals.size(); i++){
        triplets.push_back(Eigen::Triplet<double>(i,i,diagonals[i]));
    }
    mat->setFromTriplets(triplets.begin(), triplets.end());
}

void generateBxMatrix(long lattice_size, Eigen::SparseMatrix<double> *mat){
    long N = lattice_size * lattice_size;
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    std::vector<Eigen::Triplet<double>> triplets;
    for(long i = 1; i < matrixDim ;i *= 2){
        long currentPowerOfTwo = i;
        bool oneMode = true;
        long subCount = 0;
        for(long currentRow = 0; currentRow < matrixDim && currentRow + currentPowerOfTwo < matrixDim; currentRow++){
            if(oneMode){
                triplets.push_back(Eigen::Triplet<double>(currentRow, currentRow + currentPowerOfTwo, 1.0));
                triplets.push_back(Eigen::Triplet<double>(currentRow + currentPowerOfTwo, currentRow, 1.0));
            }
            subCount++;
            if(subCount == currentPowerOfTwo ){
                oneMode = !oneMode;
                subCount = 0;
            }
        }
    }
    mat->setFromTriplets(triplets.begin(), triplets.end());
}


DynamicMatrix generateAdjacentInteractionZ(long i, long j){
    DynamicMatrix output(1, 1);

    DynamicMatrix SigmaZI(2,2);
    SigmaZI.set(0, 0, 1.0);
    SigmaZI.set(1, 1, -1.0);

    DynamicMatrix I2(2,2);
    I2.set(0, 0, 1.0);
    I2.set(1, 1, 1.0);

    output.set(0, 0, 1.0);
    for (long x = 1; x < N+1; x++){
        if (x == i || x == j)
            output.tensorInPlace(SigmaZI);
        else
            output.tensorInPlace(I2);
    }
    return output;
}

void generateJMatrix(int lattice_size, Eigen::SparseMatrix<double> *mat ){
    long N = lattice_size * lattice_size;
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    JTerms = DynamicMatrix(matrixDim, matrixDim);
    for (long q = 1; q < N+1; q++){
        if (q + 1 < (N+1) && (q + 2) % lattice_size != 0)          
            JTerms.addInPlace(generateAdjacentInteractionZ(q, q+1));
        if (q + lattice_size < (N+1))
            JTerms.addInPlace(generateAdjacentInteractionZ(q, q + lattice_size));
    }
    JTerms.setEigenSparseMatrix(mat);

}


int main(int argc, *char[] argc){
    //Initialize arguments

    if(argc < 5){
        std::cout <<  "This program requires at least 5 arguments" << std::endl;
        std::cout << "Usage: ./matgen <Lattice Size in 1D> <J multiplier> <Bx Multiplier> <Bz Multiplier>" << std::endl;
        return 0;
    }

    bool verbose = false;
    int lattice_size;
    double J;
    double Bx;
    double Bz;

    try{
        std::string sizeStr(argv[1]);
        lattice_size = std::stoi(sizeStr);
    }
    catch(const std::invalid_argument & e){
        printf("Invalid argument %s supplied for lattice_size.\n", argv[1]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        printf( "Out of range argument %s supplied for lattice_size.\n" , argv[1]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string JString(argv[2]);
        J = std::stod(JString);
    }
    catch(const std::invalid_argument & e){
        printf( "Invalid argument %s supplied for J.\n", argv[2]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        printf( "Out of range argument %s supplied for J.\n" , argv[2]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BxString(argv[3]);
        Bx = std::stod(BxString);
    }
    catch(const std::invalid_argument & e){
        printf( "Invalid argument %s supplied for Bx.", argv[3]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        printf( "Out of range argument %s supplied for Bx.\n" , argv[3]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BzString(argv[4]);
        Bz = std::stod(BzString);
    }
    catch(const std::invalid_argument & e){
        printf( "Invalid argument %s supplied for Bz.\n", argv[4]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        printf( "Out of range argument %s supplied for Bz.\n" , argv[4]);
        SlepcFinalize();
        return 0;
    }

    long N = lattice_size * lattice_size;
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }

    SpMat J(matrixDim, matrixDim);
    SpMat Bx(matrixDim, matrixDim);
    Spmat Bz(matrixDim, matrixDim);

    generateJMatrix(lattice_size, &J);
    generateBxMatrix(lattice_size, &Bx);
    generateBzMatrix(lattice_size, &Bz);

    







}