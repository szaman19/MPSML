#include "IsingHamiltonian.h"

IsingHamiltonian::IsingHamiltonian(int latice_size){
    latice_size_one_dimension = latice_size;
    N = latice_size * latice_size;
    initializePauliMatrices();
    generateJMatrix();
    generateBxMatrix();
    generateBzMatrix();
    
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

DynamicMatrix IsingHamiltonian::generateSelfInteraction(char spin, int i){
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);
    for (int x = 1; x < N+1; x++)
    {
        if (x == i)
        {
            if (spin == 'x')
                output = output.tensor(SigmaXI);
            else
                output = output.tensor(SigmaZI);
        }
        else
            output = output.tensor(I2);

    }
    
    return output;
}

DynamicMatrix IsingHamiltonian::generateAdjacentInteractionZ(int i, int j){
    DynamicMatrix output(1, 1);
    output.set(0, 0, 1.0);
    for (int x = 1; x < N+1; x++)
    {
        if (x == i)
            output = output.tensor(SigmaZI);
        else if(x == j)
            output = output.tensor(SigmaZI);
        else
            output = output.tensor(I2);
    }
    return output;
}

void IsingHamiltonian::generateJMatrix(){
    int matrixDim = 1;
    for(int i = 0; i < N; i++){
        matrixDim = 2 * matrixDim;
    }
    
    JTerms = DynamicMatrix(matrixDim, matrixDim);
    for (int q = 1; q < N+1; q++)
    {
        if (q + 1 < (N+1) && (q + 2) % latice_size_one_dimension != 0)          
            JTerms = JTerms + generateAdjacentInteractionZ(q, q+1);
        if (q + latice_size_one_dimension < (N+1))
            JTerms = JTerms + generateAdjacentInteractionZ(q, q + latice_size_one_dimension);
    }
}

void IsingHamiltonian::generateBxMatrix(){
    BxTerms = DynamicMatrix(1, 1);
    BxTerms.set(0, 0, 1.0);
    for (int q = 1; q < N+1; q++)
    {
        if (q == 1)
            BxTerms = generateSelfInteraction('x', q);
        else
        {
            BxTerms = (BxTerms + generateSelfInteraction('x', q));
        }
    }

}

void IsingHamiltonian::generateBzMatrix(){
    BzTerms = DynamicMatrix(1, 1);
    BzTerms.set(0, 0, 1.0);
    for (int q = 1; q < N+1; q++)
    {
        if (q == 1)
            BzTerms = generateSelfInteraction('z', q);
        else
            BzTerms = BzTerms + generateSelfInteraction('z', q);
    }
}

DynamicMatrix IsingHamiltonian::getHamiltonian(double J, double Bx, double Bz){
    DynamicMatrix Jmul = JTerms * J;

    DynamicMatrix Bxmul = BxTerms * (Bx * -1.0);

    DynamicMatrix Bzmul = BzTerms * (Bz * -1.0);

    return Jmul + Bxmul + Bzmul;
}