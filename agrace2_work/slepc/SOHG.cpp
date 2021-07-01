#include "SlepcIsingHamiltonian.cpp"

/* SLEPc-based Object-oriented Hamiltonian Generator */ 

static char help[] = "SLEPc-based Object-oriented Hamiltonian Generator/solver (SOHG), by Andrew Grace.";

int main(int argc, char *argv[]){

    

    SlepcInitialize( &argc, &argv, (char*) 0, help);
    

    /*  Before we bother with SLEPc, we should get arguments
        Argument 1: Lattice Size
        Argument 2: J
        Argument 3: Bx
        Argument 4: Bz
    */

    if(argc < 5 || argc > 6){
        PetscPrintf(MPI_COMM_WORLD, "This program requires 5 or 6 arguments.\n");
        PetscPrintf(MPI_COMM_WORLD, "Usage: ./sohg <Lattice Size in 1D> <J multiplier> <Bx Multiplier> <Bz Multiplier> <--verbose>\n");
        PetscPrintf(MPI_COMM_WORLD, "Please try again.\n");
        SlepcFinalize();
        return 0;
    }

    /* They passed the first test, so let's get the arguments */

    bool verbose = false;
    int lattice_size;
    double J;
    double Bx;
    double Bz;

    if(argc == 6){ 
        if(!strcmp(argv[5], "--verbose"))
        verbose = true;
        else std::cout << argv[5] << std::endl;
    }

    try{
        std::string sizeStr(argv[1]);
        lattice_size = std::stoi(sizeStr);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for lattice_size.\n", argv[1]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for lattice_size.\n" , argv[1]);
        SlepcFinalize();
        return 0;
    }


    try{
        std::string JString(argv[2]);
        J = std::stod(JString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for J.\n", argv[2]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for J.\n" , argv[2]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BxString(argv[3]);
        Bx = std::stod(BxString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bx.", argv[3]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bx.\n" , argv[3]);
        SlepcFinalize();
        return 0;
    }

    try{
        std::string BzString(argv[4]);
        Bz = std::stod(BzString);
    }
    catch(const std::invalid_argument & e){
        PetscPrintf(MPI_COMM_WORLD, "Invalid argument %s supplied for Bz.\n", argv[4]);
        SlepcFinalize();
        return 0;
    }
    catch(const std::out_of_range & e){
        PetscPrintf(MPI_COMM_WORLD, "Out of range argument %s supplied for Bz.\n" , argv[4]);
        SlepcFinalize();
        return 0;
    }

    if(verbose) PetscPrintf(MPI_COMM_WORLD, "\n\033[0;32mSLEPC-based Object-oriented Hamiltonian Generator/solver.\n2021, Andrew Grace, Binghamton University.\n\033[0mInitializing SLEPC and PETSC...\n");
    if(verbose) PetscPrintf(MPI_COMM_WORLD, "\033[0;32mDone!\n\033[0m");

    {

        /* If the user has good input, start work */

        SlepcIsingHamiltonian h (lattice_size, verbose);
        h.getSolution(J, Bx, Bz);        
        
    }
    
    SlepcFinalize();



}