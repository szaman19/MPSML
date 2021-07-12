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
 
