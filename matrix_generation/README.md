# Matrix Generation
This program creates the training data for our neural networks by solving an Ising Hamiltonian for its lowest Eigenvalue and Eigenvector. 

## How to use
- Make the program. You may need to adjust the makefile to suit your configuration. The one in this directory is currently openHPC's SLEPC and PETSC configuration.
- Run the program. It's called dhogs, or Distributed Hamiltonian Object-oriented Generator and Solver.
- The current version takes these arguments:
    - Lattice size: the one-dimensional size of the matrix (this program supports only square lattice structures).
    - J value: the value by which the J component of the Hamiltonian should be multiplied.
    - Bx value: the value by which the Bx component of the Hamiltonian should be multiplied.
    - Bz value: the value by which the Bz component of the Hamiltonian should be multiplied.
    - --verbose: this option will force the program to provide more details on the solve.
    - --forcegen: this option forces the matrix to ignore current files and regenerate all matrix components

- The next version will support processing different combinations of J, Bx, and Bz simultanenously, and it will take these arguments:
   - Lattice size: same as previous
   - CSV filename: The CSV which contains the combinations of J, Bx, and Bz.
   - --verbose: same as previous
 
 ## CSV Format
 This program will expect it's data in this format:
 J1, Bx1, Bz1
 J2, Bx2, Bz2
 .., .. , ..

 Where these values are all ASCII-printed doubles.
 
 ## Output Eigenpair format
 The output file format of this program stores one eigenpair per file. Here is the file format:
 - File format version: Little Endian 1. 4 bytes
 - Big Endian Indicator: stored as little endian regardless of system. 1 for big, 0 for little. 4 bytes
 - Eigenvalue: double, endianness depends on indicator. 8 bytes
 - Number of values in eigenvector: integer, endianness depends on indicator, 4 bytes
 - Eigenvector values. Double, endianness depends on indicator, 8 bytes each.

We currently have a complete C++ class that can load/store the eigenpair to files, and a beta python script. In the future, I will write a Python addon that uses the C++ class, as I like C++ better.

 ## How this program works.
 Due to a number of nightmares during development, my program performs a long process to generate matrices.
 - If necessary, my code will generate the J, Bx, and Bz matrices and save them to the folder that the program was compiled into. It saves them in a format that PETSC can read.
 - When done generating, my code instructs PETSC to load the matrices and perform arithmetic on them to get the matrix specified by your input J, Bx, Bz values.
 - That will run into SLEPC, which performs the eigenvector calculation. It will temporarily save the eigenvector as a native PETSC vector to the filesystem
 - My code will combine that vector with the eigenvalue, package it into the nice eigenpair format, delete the old file, and place the new one into the directory.
 
 
