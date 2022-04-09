#include "ArgumentMethods.h"
#include "IsingHamiltonian.h"
#include "StringsAndFormats.h"
#include "debug_flags.hpp"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>

static char petsc_help[] = "PetscHelp";

#define dprintf(text) PetscPrintf(MPI_COMM_WORLD, text);

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);
	SlepcInitialize(&argc, &argv, (char *)0, petsc_help);
	int num_qbits;
	int debug_arg_start = 0;

	std::vector<IsingParmStruct> parameters_to_test;

	if (argc <= 2) {
		std::cout << EXECUTABLE_NAME << ERR_INSUFFICIENT_ARGS << std::endl;
		std::cout << USAGE_RANGE << std::endl;
		std::cout << USAGE_CSV << std::endl;
		return -1;
	}

	num_qbits = get_number_from_string<int>(std::string(argv[1]), "reading number of qubits");

	if (argv[2] == std::string("csv")) {
		if (argc >= 4) {
			debug_arg_start    = 4;
			parameters_to_test = get_parms_from_csv(std::string(argv[3]));
		} else {
			std::cout << "CSV Mode" << ERR_INVALID_ARGS << std::endl;
			std::cout << USAGE_CSV << std::endl;
			return -1;
		}

	} else if (argv[2] == std::string("range")) {
		if (argc >= 10) {
			debug_arg_start    = 10;
			double j           = get_number_from_string<double>(std::string(argv[3]));
			int    num_bx      = get_number_from_string<int>(std::string(argv[4]));
			double init_bx     = get_number_from_string<double>(std::string(argv[5]));
			double stop_bx     = get_number_from_string<double>(std::string(argv[6]));
			int    num_bz      = get_number_from_string<int>(std::string(argv[7]));
			double init_bz     = get_number_from_string<double>(std::string(argv[8]));
			double stop_bz     = get_number_from_string<double>(std::string(argv[9]));
			parameters_to_test = get_parms_from_range(j, num_bx, init_bx, stop_bx, num_bz, init_bz, stop_bz);
		} else {
			std::cout << "Range Mode" << ERR_INVALID_ARGS << std::endl;
			std::cout << USAGE_RANGE << std::endl;
			return -1;
		}
	} else {
		std::cout << EXECUTABLE_NAME << ERR_INVALID_MODE << std::endl;
		std::cout << USAGE_RANGE << std::endl;
		std::cout << USAGE_CSV << std::endl;
		return -1;
	}

	debug_flags flags;
	process_debug_args(debug_arg_start, argc, argv, &flags);
	{
		IsingHamiltonian hamiltonian(num_qbits, flags);

		auto start_solve_time = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < parameters_to_test.size(); i++) {
			IsingParmStruct current_parameter = parameters_to_test[i];

			hamiltonian.solve_for_eigenvector(current_parameter.j_scalar, current_parameter.bx_scalar, current_parameter.bz_scalar);
		}

		auto end_solve_time              = std::chrono::high_resolution_clock::now();
		auto duration_for_this_solve     = std::chrono::duration_cast<std::chrono::milliseconds>(end_solve_time - start_solve_time);
		long milliseconds_for_this_solve = (long)duration_for_this_solve.count();

		if (flags.run_performance_metrics) {
			PetscPrintf(PETSC_COMM_WORLD, "Performance Metrics :\n");
			PetscPrintf(PETSC_COMM_WORLD, "J Generation Time   : %ld ms\n", hamiltonian.j_generation_time);
			PetscPrintf(PETSC_COMM_WORLD, "Bx Generation Time  : %ld ms\n", hamiltonian.bx_generation_time);
			PetscPrintf(PETSC_COMM_WORLD, "Bz Generation Time  : %ld ms\n", hamiltonian.bz_generation_time);
			PetscPrintf(PETSC_COMM_WORLD, "Average Solve Time  : %ld ms\n", milliseconds_for_this_solve / parameters_to_test.size());
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	SlepcFinalize();
	MPI_Finalize();
	return 0;
}
