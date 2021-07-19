#include "IsingHamiltonian.h"
#include <chrono>
#include<iostream>

IsingHamiltonian::IsingHamiltonian(long latice_size, int which){
    latice_size_one_dimension = latice_size;
    N = latice_size * latice_size;
    initializePauliMatrices();
    std::string JMatFileName = "JMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::string BzMatFileName = "BzMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::string BxMatFileName = "BxMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::ifstream jfile(JMatFileName);
    std::ifstream bzfile(BzMatFileName);
    std::ifstream bxfile(BxMatFileName);

    if(!jfile.good() && which == 0){
        auto start = std::chrono::system_clock::now();
        generateJMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        JTerms.savePetsc(JMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "J Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }
    if(!bzfile.good() && which == 1){
        auto start = std::chrono::system_clock::now();
        generateBzMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        BzTerms.savePetsc(BzMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "Bz Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }
    if(!bxfile.good() && which == 2){
        auto start = std::chrono::system_clock::now();
        generateBxMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        BxTerms.savePetsc(BxMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "Bx Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }

    jfile.close();
    bzfile.close();
    bxfile.close();



}

IsingHamiltonian::IsingHamiltonian(long latice_size, int which, bool sv){
    latice_size_one_dimension = latice_size;
    N = latice_size * latice_size;
    initializePauliMatrices();
    std::string JMatFileName = "JMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::string BzMatFileName = "BzMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::string BxMatFileName = "BxMat_" + std::to_string(latice_size) + ".dynamicmatrix";
    std::ifstream jfile(JMatFileName);
    std::ifstream bzfile(BzMatFileName);
    std::ifstream bxfile(BxMatFileName);

    if(!jfile.good() && which == 0){
        auto start = std::chrono::system_clock::now();
        generateJMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        if(sv)
            JTerms.savePetsc(JMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "J Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }
    if(!bzfile.good() && which == 1){
        auto start = std::chrono::system_clock::now();
        generateBzMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        if(sv)
            BzTerms.savePetsc(BzMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "Bz Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }
    if(!bxfile.good() && which == 2){
        auto start = std::chrono::system_clock::now();
        generateBxMatrixExperimental();
        auto save = std::chrono::system_clock::now();
        if(sv)
        BxTerms.savePetsc(BxMatFileName);
        auto end = std::chrono::system_clock::now();
        auto totalElapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto genElapsed = std::chrono::duration_cast<std::chrono::seconds>(save - start);
        auto saveElapsed =std::chrono::duration_cast<std::chrono::seconds>(end - save);
        std::cout << "Bx Matrix times:" << std::endl << "Generation: "  << genElapsed.count() << std::endl << "Save: " << saveElapsed.count() << std::endl << "Total: " << totalElapsed.count() << std::endl << std::endl;
    }

    jfile.close();
    bzfile.close();
    bxfile.close();


}





void IsingHamiltonian::initializePauliMatrices(){
    SigmaXI = DynamicMatrix(2,2);
    SigmaXI.set(0, 1, 1.0);
    SigmaXI.set(1, 0, 1.0);

    SigmaZI = DynamicMatrix(2,2);
    SigmaZI.set(0, 0, 1.0);
    SigmaZI.set(1, 1, -1.0);
    
    I2 = DynamicMatrix(2,2);
    I2.set(0, 0, 1.0);
    I2.set(1, 1, 1.0);
}

DynamicMatrix IsingHamiltonian::generateSelfInteraction(char spin, long i){
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);
    for (long x = 1; x < N+1; x++){
        if (x == i){
            if (spin == 'x')
                output.tensorInPlace(SigmaXI);
            else
                output.tensorInPlace(SigmaZI);
        }
        else
            output.tensorInPlace(I2);
    }
    return output;
}

void IsingHamiltonian::generateBxMatrixExperimental(){
    //Note
    //Row 0 of Bx: 1 only when column is a power of 2, including 1
    //Col 0 of Bx: 1 only when row is power of 2, including 1
    //0 along diagonal
    //From Row 0: 1 along diagonal where column is 1. There are as many ones as the column number. Then there are that many zeros
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    
    BxTerms = DynamicMatrix(matrixDim,matrixDim);

    //Calculate nonzeros required
    long nonzeros = N * matrixDim;
    BxTerms.reserve(nonzeros);


    long startPoint = 1;
    auto begin = std::chrono::high_resolution_clock::now();
    for(long i = 1; i < matrixDim ;i *= 2){
        long currentPowerOfTwo = i;
        bool oneMode = true;
        long subCount = 0;
        for(long currentRow = 0; currentRow < matrixDim && currentRow + currentPowerOfTwo < matrixDim; currentRow++){
            if(oneMode){
                BxTerms.set(currentRow, currentRow + currentPowerOfTwo, 1.0);
                BxTerms.set(currentRow + currentPowerOfTwo, currentRow, 1.0);
            }
            subCount++;
            if(subCount == currentPowerOfTwo ){
                oneMode = !oneMode;
                subCount = 0;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Bx: " << duration.count() << " seconds." << std::endl;

}

void IsingHamiltonian::generateBzMatrixExperimental(){
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    
    BzTerms = DynamicMatrix(matrixDim,matrixDim);
    
    //Levels of pattern recursion: lattice_size
    double start = N;
    std::vector<double> diagonals = generateBzPattern(1,N);


    for(int i = 0; i < diagonals.size(); i++){
        BzTerms.set(i,i,diagonals[i]);
    }

    
}


std::vector<double> IsingHamiltonian::generateBzPattern(int level, double input){

    if(level == N){
        std::vector <double> returner;
        returner.push_back(input);
        returner.push_back(input - 2);

        return returner;
    }
    else{
        std::vector<double> returner;
        
        std::vector<double> r1 = generateBzPattern(level + 1, input);
        std::vector<double> r2 = generateBzPattern(level + 1, input - 2);
        
        returner.insert(std::end(returner), std::begin(r1), std::end(r1));
        returner.insert(std::end(returner), std::begin(r2), std::end(r2));

        return returner;
    }


}





DynamicMatrix IsingHamiltonian::generateAdjacentInteractionZ(long i, long j){
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);
    for (long x = 1; x < N+1; x++){
        if (x == i || x == j)
            output.tensorInPlace(SigmaZI);
        else
            output.tensorInPlace(I2);
    }
    return output;
}



void IsingHamiltonian::generateAdjacentBeta(long matdim, long i, long j, std::vector<double> &toAdd){
    std::vector<double> output;
    for(long i = 0; i < matdim ; i++){
        output.push_back(0.0);
    }
    //q and q+1
    output[0] = 1.0;

    for (long x = 1; x < N+1; x++){
        long twoRaisedX = 1;
        for(int i = 0; i < x; i++){
            twoRaisedX *=2;
        }

        long thisStep = matdim / (twoRaisedX / 2);
        long nextStep = matdim / twoRaisedX;
    
        if (x == i || x == j){
            
            for(long i = 0; i < output.size(); i+= thisStep){
                output[i + nextStep] = -1.0 * output[i];
            }
           
        }
        else{
            for(long i = 0; i < output.size(); i+= thisStep){
                output[i + nextStep] = output[i];
            }
        }  
    }

    bool writeMode = false;
    if(toAdd.size() < output.size()) writeMode = true;
    for(long i = 0; i < output.size(); i++){
        if(writeMode) toAdd.push_back(output[i]);
        else toAdd.at(i) += output[i];
    }



}

void IsingHamiltonian::generateJMatrixExperimental(){
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }

    std::vector<double> diagonal;
    JTerms = DynamicMatrix(matrixDim, matrixDim);

    for (long q = 1; q < N+1; q++){
        std::vector<double> addedDiagonal;
        
        if (q + 1 < (N+1) && (q + 2) % latice_size_one_dimension != 0)
            generateAdjacentBeta(matrixDim, q, q+1, addedDiagonal);
        if (q + latice_size_one_dimension < (N+1))
            generateAdjacentBeta(matrixDim, q, q + latice_size_one_dimension, addedDiagonal);

        bool writeMode = addedDiagonal.size() > diagonal.size();
        for(long i = 0; i < addedDiagonal.size(); i++){
            if(writeMode) diagonal.push_back(addedDiagonal[i]);
            else diagonal[i] += addedDiagonal[i];
        }
    }

    //Write the diagonal to JTerms
    for(long i = 0; i < matrixDim && i < diagonal.size(); i++){
        JTerms.set(i,i,diagonal[i]);
    }


}

void IsingHamiltonian::generateJMatrix(){
    long matrixDim = 1;
    for(long i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    
    JTerms = DynamicMatrix(matrixDim, matrixDim);
    for (long q = 1; q < N+1; q++){
        if (q + 1 < (N+1) && (q + 2) % latice_size_one_dimension != 0)          
            JTerms.addInPlace(generateAdjacentInteractionZ(q, q+1));
        if (q + latice_size_one_dimension < (N+1))
            JTerms.addInPlace(generateAdjacentInteractionZ(q, q + latice_size_one_dimension));
    }
}

void IsingHamiltonian::generateBxMatrix(){
    BxTerms = DynamicMatrix(1, 1);
    BxTerms.set(0, 0, 1.0);
    for (long q = 1; q < N+1; q++){
        if (q == 1)
            BxTerms = generateSelfInteraction('x', q);
        else
            BxTerms.addInPlace(generateSelfInteraction('x', q));
    }

}

void IsingHamiltonian::generateBzMatrix(){
    BzTerms = DynamicMatrix(1, 1);
    BzTerms.set(0, 0, 1.0);
    for (long q = 1; q < N+1; q++){
        if (q == 1)
            BzTerms = generateSelfInteraction('z', q);
        else
            BzTerms.addInPlace(generateSelfInteraction('z', q));
    }
}

DynamicMatrix IsingHamiltonian::getHamiltonian(double J, double Bx, double Bz){
    DynamicMatrix Jmul = JTerms * J;
    DynamicMatrix Bxmul = BxTerms * (Bx * -1.0);
    DynamicMatrix Bzmul = BzTerms * (Bz * -1.0);
    return Jmul + Bxmul + Bzmul;
}