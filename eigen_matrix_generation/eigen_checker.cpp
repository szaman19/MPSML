#include<iostream>
#include<Eigen/Dense>
#include<Eigen/KroneckerProduct>

using namespace Eigen;

Matrix2d I2;
Matrix2d SigmaX;
Matrix2d SigmaZ;
long lattice_size;

MatrixXd generateSelfInteraction(char spin, long i){
    MatrixXd output(1,1);
    output << 1;
    for(long x = 1; x < (lattice_size * lattice_size) + 1; x++){
        if(x == i){
            if(spin == 'x'){
                output = kroneckerProduct(output, SigmaX).eval();
            }
            else{
                output = kroneckerProduct(output, SigmaZ).eval();
            }
        }
        else{
            output = kroneckerProduct(output, I2).eval();
        }

    }
    return output;
}

MatrixXd generateAdjacentInteraction(long i, long j){
    MatrixXd output(1,1);
    output << 1;

    for(long x = 1; x < (lattice_size * lattice_size) + 1; x++){
        if(x == i || x == j){
            output = kroneckerProduct(output, SigmaZ).eval();
        }
        else{
            output = kroneckerProduct(output, I2).eval();
        }
    }
    return output;
}

MatrixXd genJMat(){
    long matrixDim = 1;
    for(long i = 0; i < (lattice_size * lattice_size); i++){
        matrixDim = 2 * matrixDim;
    }

    MatrixXd JTerms(matrixDim, matrixDim);
    long N = (lattice_size * lattice_size);
    for (long q = 1; q < (lattice_size * lattice_size)+1; q++){
        if (q + 1 < (N+1) && (q + 2) % lattice_size != 0)          
            JTerms += generateAdjacentInteraction(q, q+1);
        if (q + lattice_size < (N+1))
            JTerms += generateAdjacentInteraction(q, q + lattice_size);
    }

    return JTerms;
}

MatrixXd genBzMat(){
    long matrixDim = 1;
    for(long i = 0; i < (lattice_size * lattice_size); i++){
        matrixDim = 2 * matrixDim;
    }

    MatrixXd BzTerms(matrixDim, matrixDim);
    for (long q = 1; q < (lattice_size * lattice_size)+1; q++){
        if (q == 1)
            BzTerms = generateSelfInteraction('z', q);
        else
            BzTerms += (generateSelfInteraction('z', q));
    }
    return BzTerms;

}

MatrixXd genBxMat(){
    long matrixDim = 1;
    for(long i = 0; i < (lattice_size * lattice_size); i++){
        matrixDim = 2 * matrixDim;
    }

    MatrixXd BxTerms(matrixDim, matrixDim);
    for (long q = 1; q < (lattice_size * lattice_size)+1; q++){
        if (q == 1)
            BxTerms = generateSelfInteraction('x', q);
        else
            BxTerms += (generateSelfInteraction('x', q));
    }
    return BxTerms;

}




int main(){

    double J, Bx, Bz;

    lattice_size = 2;
    J = 1.0;
    Bx = 1.0;
    Bz = 1.0;
    

    I2 << 1, 0,
          0, 1;
    
    SigmaX << 0, 1,
              1, 0;

    SigmaZ << 1, 0,
              0, -1;

    MatrixXd JMat = genJMat();
    MatrixXd BzMat = genBzMat();
    MatrixXd BxMat = genBxMat();

    std::cout << JMat << std::endl <<std::endl;
    std::cout << BzMat <<std::endl <<std::endl;
    std::cout << BxMat << std::endl <<std::endl;





}
