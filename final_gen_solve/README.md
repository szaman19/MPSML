# How to use

## Setting everything up.
- Run ``make dirs`` to make the required directories for this program.
- Run ``make gencsv `` to build gencsv, a CSV generator that makes input files for the generator.
- Generate the CSV with ``./gencsv <J value> <Stops in Bx> <Stops in Bz> <Initial Bx>  <Initial Bz> <Stop Bx> <Stop Bz> <file name>``
- Run ``make`` to make the generator
- Modify the SLURM script to run the file with your parameters. The final line in that file should look like ``srun --mpi=pmix_v2 ./matgen <Lattice size in 1D> <CSV file>  <PETSC args>``
- The eigenvectors will be stored in a single file called results.eigenset. You may use the Python library or the C++ library to read these values
 
 ## Python Eigenset Reader Documentation
 The Python Eigenset reader lets you access the eigenvalues in the .eigenset file. Each eigenvector/eigenvalue is stored in a IsingEigenset object, which contains the following values:
 - J
 - Bx
 - Bz
 - Eigenvalue
 - Eigenvector[], a list of doubles that represents the eigenvector.

The Eigenset object contains a list of these IsingEigensets named eigenpairs. It also contains methods to read eigenset files. Here's a demo program that prints out the J, Bx, and Bz files for a .eigenset file:

    import EigensetReader
    //Create the Eigenset and read the file
    e = EigensetReader.Eigenset()
    e.read("results.eigenset")
    
    for x in range(0, e.numberEigenvectors):
        print(str(e.eigenpairs[x].J) + " " + str(e.eigenpairs[x].Bx) + " " + str(e.eigenpairs[x].Bz) + " " + str(e.eigenpairs[x].Eigenvalue) )
        
