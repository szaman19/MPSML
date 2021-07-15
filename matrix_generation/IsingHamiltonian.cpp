#include "IsingHamiltonian.h"

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
        generateJMatrix();
        JTerms.savePetsc(JMatFileName);
        //JTerms.save(JMatFileName + "e");
    }
    if(!bzfile.good() && which == 1){
        generateBzMatrixExperimental();
        BzTerms.savePetsc(BzMatFileName);
        //BzTerms.save(BzMatFileName + "e");
    }
    if(!bxfile.good() && which == 2){
        generateBxMatrixExperimental();
        BxTerms.savePetsc(BxMatFileName);
        //std::cout << BxTerms.printLatex();
        //BxTerms.save(BxMatFileName + "e");


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
    std::vector<long> powersOfTwo;
    long startPoint = 1;
    while(matrixDim % startPoint == 0){
        powersOfTwo.push_back(startPoint);
        startPoint *= 2;
    }

    //Generate the upper half of the matrix
    for(long i = 0; i < powersOfTwo.size();i++){
        long currentPowerOfTwo = powersOfTwo[i];
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
}

void IsingHamiltonian::generateBzMatrixExperimental(){
    BzTerms = DynamicMatrix(1,1);

    //Note
    //We're only going to hit tensorInPlace(SigmaZI) when i == x

    for(long i = 1; i < N+1; i++ ){
        DynamicMatrix temp(1,1);
        temp.set(0,0,1.0);
        //generateSelfInteraction(z,i)
        for(long x = 1; x < N+1; x++){
            if(x == i) temp.tensorInPlace(SigmaZI);
            else temp.tensorInPlace(I2);
        }
        if(i == 1) BzTerms = temp;
        else BzTerms.addInPlaceRef(temp);
        BzTerms.clean();
        
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