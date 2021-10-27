#include<iostream>
#include<Eigen/Dense>
#include<Eigen/KroneckerProduct>
#include"DynamicMatrix.h"
#include<cmath>

using namespace Eigen;
class EigenIsingHamiltonian{
    public:

    Matrix2d I2;
    Matrix2d SigmaX;
    Matrix2d SigmaZ;

    MatrixXd Jmat;
    MatrixXd Bzmat;
    MatrixXd Bxmat;

    long lattice_size;  

    EigenIsingHamiltonian(long lattice_size){
        this->lattice_size = lattice_size;
        I2 << 1, 0,
              0, 1;
    
        SigmaX << 0, 1,
                  1, 0;

        SigmaZ << 1, 0,
                  0,-1;

        Jmat = genJMat();
        Bxmat = genBxMat();
        Bzmat = genBzMat();

        
    }


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

    bool compareToDynamicMatrix(DynamicMatrix &other, int whichMatrix){
        MatrixXd * pointer;
        double epsilon = 0.000001;
        if(whichMatrix == 0) pointer = &Jmat;
        else if(whichMatrix == 1) pointer = &Bzmat;
        else pointer = &Bxmat;

        bool same = true;
        for(int i = 0; i < pointer->rows() && same; i++){
            for(int j = 0; j < pointer->cols() && same; j++){
                if(pointer->operator()(i,j) == other.get(i,j)){

                }
                else if(std::isnan(pointer->operator()(i,j))){


                }
                else if(std::abs(pointer->operator()(i,j) - other.get(i,j)) < epsilon ){

                }

                else same = false;
            }
        }
        return same;

    }






};