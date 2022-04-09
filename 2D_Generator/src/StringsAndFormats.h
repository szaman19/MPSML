#ifndef PHS_STRINGS
#define PHS_STRINGS 1

const std::string EXECUTABLE_NAME = "matgen";
const std::string USAGE_RANGE = "Usage: ./" + EXECUTABLE_NAME + " <num_qbits> range <j> <num_bx> <init_bx> <stop_bx> <num_bz> <init_bz> <stop_bz>";
const std::string USAGE_CSV = "Usage: ./" + EXECUTABLE_NAME + " <num_qbits> csv <filename>";
const std::string DEBUG_ARGS = "Available debug flags:\n"
"--save-components-ascii        : saves J, Bx, Bz as human-readable files.\n"
"--save-components-petsc        : saves J, Bx, Bz as PETSC-readable files.\n"
"--print-eigenvalue-information : Prints additional info for each eigenvalue.\n"
"--run-performance-metrics     : show performance metrics after solve.\n"
"--perform-eigenvalue-check    : ensure the matrix was diagonalized correctly by checking that mat * eigvec = lambda * eigvec.\n"
"--verbose                     : print information during matrix generation.\n"
"--file=<filestr>              : saves eigenset output to a different file.\n";

const std::string ERR_INSUFFICIENT_ARGS = " requires more arguments to solve the Ising Hamiltonian.";
const std::string ERR_INVALID_MODE = " does not currently support this mode";
const std::string ERR_INVALID_ARGS = " was provided either invalid or an invalid number of arguments";

#endif