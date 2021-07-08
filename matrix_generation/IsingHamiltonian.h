#pragma once
#include "DynamicMatrix.h"
#include<iostream>
//  ============================================================
//  Ising Hamiltonian Generator
//  Created by Andrew T Grace, 2021
//  Binghamton University, Computer Science Research
//  Overseen by Professors Kenneth Chiu, Wei-Cheng Lee.
//
//  This code creates an Ising hamiltonian, given your specified
//  J, Bx, and Bz scalar terms. It returns it as a DynamicMatrix
//  which allows us to return matrices in different formats for
//  different eigensolvers. We can use this to train our neural 
//  network for solving Ising hamiltonians with more particles.
//  ============================================================ 



class IsingHamiltonian{
    
    private:

        //Matrix storage, generated during initialization
        

        //Pauli Matrix Storage
        DynamicMatrix SigmaXI;
        DynamicMatrix SigmaZI;
        DynamicMatrix I2;

        //Functions for getting Pauli matrices, I2 = [1,0;0,1]
        void initializePauliMatrices();

        //Helper functions for finding interactions
        DynamicMatrix generateSelfInteraction(char, long);
        DynamicMatrix generateAdjacentInteractionZ(long, long);

        //Functions for initializing JTerms, BxTrems
        void generateJMatrix();
        void generateBxMatrix();
        void generateBzMatrix();

        //N should always equal latice_size_one_dimension^2
        long N;
        long latice_size_one_dimension;

    public:
        IsingHamiltonian(long, int);
        DynamicMatrix getHamiltonian(double, double, double);

        DynamicMatrix JTerms;
        DynamicMatrix BxTerms;
        DynamicMatrix BzTerms;


        
    
        

};