# How to use

## Setting everything up.
- Run ``make dirs`` to make the required directories for this program.
- Run ``make gencsv `` to build gencsv, a CSV generator that makes input files for the generator.
- Generate the CSV with ``./gencsv <J value> <Stops in Bx> <Stops in Bz> <Initial Bx> <Stop Bx> <Initial Bz> <Stop Bz> <file name>``
- Run ``make`` to make the generator
- Modify the SLURM script to run the file with your parameters. The final line in that file should look like ``srun --mpi=pmix_v2 ./matgen <Lattice size in 1D> file <CSV file> <--verbose, --validate, or --force-gen> <PETSC args>``
- After executing, the eigenvector files should be stored in ``./generated_eigenvectors``. You can read the file by using ``vectorlib.cpp`` in your programs to read the file into a PETSCVectorLoader class.
  - Make an instance of a PETSCVectorLoader
  - Use its ``read(filename)`` function to read the eigenvector file into the fields.
  - Use ``std::vector<double>& getVals() `` and ``double getEigenvalue`` to read the eigenvalues.
