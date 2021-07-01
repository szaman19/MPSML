#include "DynamicMatrix.h"
#include "IsingHamiltonian.h"
#include<iostream>

//  ============================================================
//  Ising Hamiltonian Generator Driver
//  Created by Andrew T Grace, 2021
//  Binghamton University, Computer Science Research
//  Overseen by Professors Kenneth Chiu, Wei-Cheng Lee.
//
//  This code uses the DynamicMatrix and IsingHamiltoninan 
//  classes to test out my code.
//  ============================================================

int main(){
    IsingHamiltonian IH(2);
    std::cout << "J Matrix:" << std::endl << IH.getHamiltonian(1.0, 0.0, 0.0) << std::endl;
    std::cout << "Bx Matrix:" << std::endl << IH.getHamiltonian(0.0, 1.0, 0.0) << std::endl;
    std::cout << "Bz Matrix:" << std::endl << IH.getHamiltonian(0.0, 0.0, 1.0) << std::endl;
    std::cout << "Full Hamiltonian (all terms 1):" << std::endl << IH.getHamiltonian(1.0, 1.0, 1.0) << std::endl;
    return 0;
}