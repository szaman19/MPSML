#include "IsingHamiltonian.h"

IsingHamiltonian::IsingHamiltonian(int num_qbits) {
	this->num_qbits          = num_qbits;
	this->num_diagonal_terms = 1 << num_qbits;
	this->num_diagonal_qbits = (int)sqrt(num_qbits);
	generate_bx();
	generate_bz();
	generate_j();
}

IsingHamiltonian::IsingHamiltonian(int num_qbits, debug_flags flags) {
	this->num_qbits          = num_qbits;
	this->num_diagonal_terms = 1 << num_qbits;
	this->num_diagonal_qbits = (int)sqrt(num_qbits);
	this->flags              = flags;

	auto start_generation = std::chrono::high_resolution_clock::now();
	vprint("Generating Bx\n");
	generate_bx();
	auto end_bx_generation = std::chrono::high_resolution_clock::now();
	vprint("Generating Bz\n");
	generate_bz();
	auto end_bz_generation = std::chrono::high_resolution_clock::now();
	vprint("Generating J\n");
	if(flags.generate_1d){
		generate_j_1d();
	}
	else{
		generate_j();
	}
	
	vprint("Finished matrix generation.\n") auto end_j_generation = std::chrono::high_resolution_clock::now();
	bx_generation_time                                            = std::chrono::duration_cast<std::chrono::milliseconds>(end_bx_generation - start_generation).count();
	bz_generation_time                                            = std::chrono::duration_cast<std::chrono::milliseconds>(end_bz_generation - end_bx_generation).count();
	j_generation_time                                             = std::chrono::duration_cast<std::chrono::milliseconds>(end_j_generation - end_bz_generation).count();

	if (flags.save_component_matrices_ascii == true) {
		std::string suffix = std::to_string(num_qbits) + "qbits.petscmat";
		save_matrix(&Bx, "Bx_" + suffix, PETSC_ASCII);
		save_matrix(&Bz, "Bz_" + suffix, PETSC_ASCII);
		save_matrix(&J, "J_" + suffix, PETSC_ASCII);
	}
	if (flags.save_component_matrices_petsc == true) {
		std::string suffix = std::to_string(num_qbits) + "qbits.txt";
		save_matrix(&Bx, "Bx_" + suffix, PETSC_PROPRIETARY);
		save_matrix(&Bz, "Bz_" + suffix, PETSC_PROPRIETARY);
		save_matrix(&J, "J_" + suffix, PETSC_PROPRIETARY);
	}
}

IsingHamiltonian::~IsingHamiltonian() {
	PW(MatDestroy(&Bx));
	PW(MatDestroy(&Bz));
	PW(MatDestroy(&J));
}

void IsingHamiltonian::solve_for_eigenvector(double j_scalar, double bx_scalar, double bz_scalar) {
	// auto start_solve_time = std::chrono::high_resolution_clock::now();
	Mat Sum;
	PW(MatCreate(PETSC_COMM_WORLD, &Sum));
	PW(MatSetSizes(Sum, PETSC_DECIDE, PETSC_DECIDE, num_diagonal_terms, num_diagonal_terms));
	PW(MatSetFromOptions(Sum));
	PW(MatSetUp(Sum));
	PW(MatZeroEntries(Sum));

	PW(MatAssemblyBegin(Sum, MAT_FINAL_ASSEMBLY));
	PW(MatAssemblyEnd(Sum, MAT_FINAL_ASSEMBLY));
	vprint("Beginning solve for j=%lf, bx=%lf, bz=%lf\n", j_scalar, bx_scalar, bz_scalar);
	PW(MatAXPY(Sum, -1.0 * j_scalar, J, DIFFERENT_NONZERO_PATTERN));
	PW(MatAXPY(Sum, -1.0 * bx_scalar, Bx, DIFFERENT_NONZERO_PATTERN));
	PW(MatAXPY(Sum, -1.0 * bz_scalar, Bz, DIFFERENT_NONZERO_PATTERN));
	/* Set up SLEPC */
	Vec         imaginary_eigenvector, real_eigenvector;
	PetscScalar imaginary_eigenvalue, real_eigenvalue;
	EPS         solver;
	int         nconv, size;

	/* Initialize vectors */
	PW(MatCreateVecs(Sum, NULL, &imaginary_eigenvector));
	PW(MatCreateVecs(Sum, NULL, &real_eigenvector));
	PW(EPSCreate(PETSC_COMM_WORLD, &solver));
	PW(EPSSetOperators(solver, Sum, NULL));
	PW(EPSSetProblemType(solver, EPS_HEP));
	PW(EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL));
	PW(EPSSetFromOptions(solver));
	PW(EPSSolve(solver));
	PW(EPSGetConverged(solver, &nconv));
	PW(EPSGetEigenpair(solver, 0, &real_eigenvalue, &imaginary_eigenvalue, real_eigenvector, imaginary_eigenvector));
	PW(EPSDestroy(&solver));

	// If the user requested eigenvalue verification, do it here.
	if (flags.perform_eigenvalue_check) {
		Vec actual_multiplied_eigenvector;
		Vec expected_multiplied_eigenvector;
		PW(MatCreateVecs(Sum, NULL, &actual_multiplied_eigenvector));
		PW(MatCreateVecs(Sum, NULL, &expected_multiplied_eigenvector));
		PW(MatMult(Sum, real_eigenvector, actual_multiplied_eigenvector));
		PW(VecCopy(real_eigenvector, expected_multiplied_eigenvector));
		PW(VecScale(expected_multiplied_eigenvector, real_eigenvalue));
		PW(VecAXPY(actual_multiplied_eigenvector, -1.0, expected_multiplied_eigenvector));
		PW(VecPointwiseMult(actual_multiplied_eigenvector, actual_multiplied_eigenvector, actual_multiplied_eigenvector));
		PetscScalar average_squared_error;
		PW(VecSum(actual_multiplied_eigenvector, &average_squared_error));
		average_squared_error = average_squared_error / size;

		if (average_squared_error > flags.eigenvalue_tolerance) {
			PW(PetscPrintf(PETSC_COMM_WORLD, "Eigenvalue for j=%f, bx=%f, bz=%f failed eigenvector test. Eigenvalue = %f, error = %f.\n", j_scalar, bx_scalar, bz_scalar, real_eigenvalue, average_squared_error));
		}
		PW(VecDestroy(&actual_multiplied_eigenvector));
		PW(VecDestroy(&expected_multiplied_eigenvector));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	sync_save_eigenvector(flags.filename, &real_eigenvector, real_eigenvalue, j_scalar, bx_scalar, bz_scalar);
	MPI_Barrier(MPI_COMM_WORLD);

	PW(VecDestroy(&imaginary_eigenvector));
	PW(VecDestroy(&real_eigenvector));
	PW(MatDestroy(&Sum));

	MPI_Barrier(MPI_COMM_WORLD);
}

// Generate adjacent interaction terms using a faster, deconstructed version of the tensor product with the Z .
void IsingHamiltonian::generate_adjacent_interaction_petsc(PetscInt i, PetscInt j, Vec *diagonal_vector) {
	i -= 1;
	j -= 1;
	PetscInt     local_processor_vector_begin;
	PetscInt     local_processor_vector_end;
	PetscScalar *local_processor_array_pointer;
	PW(VecGetOwnershipRange(*diagonal_vector, &local_processor_vector_begin, &local_processor_vector_end));
	PW(VecGetArray(*diagonal_vector, &local_processor_array_pointer));

	for (PetscInt x = local_processor_vector_begin; x < local_processor_vector_end; x++) {
		PetscInt    i_if_id = ((x % (1 << (num_qbits - i))) < (1 << (num_qbits - i - 1))) ? 1 : -1;
		PetscInt    j_if_id = ((x % (1 << (num_qbits - j))) < (1 << (num_qbits - j - 1))) ? 1 : -1;
		PetscScalar here    = (i_if_id ^ j_if_id) ? -1.0 : 1.0;
		local_processor_array_pointer[x - local_processor_vector_begin] += here;
	}

	PW(VecRestoreArray(*diagonal_vector, &local_processor_array_pointer));
	MPI_Barrier(MPI_COMM_WORLD);
}

void IsingHamiltonian::generate_j() {
	PW(MatCreate(PETSC_COMM_WORLD, &J));
	PW(MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, num_diagonal_terms, num_diagonal_terms));
	PW(MatSetFromOptions(J));
	PW(MatMPIAIJSetPreallocation(J, 1, NULL, 0, NULL));

	Vec diagonal;
	PW(VecCreate(PETSC_COMM_WORLD, &diagonal));
	PW(VecSetSizes(diagonal, PETSC_DECIDE, num_diagonal_terms));
	PW(VecSetFromOptions(diagonal));
	PW(VecZeroEntries(diagonal));
	PW(VecAssemblyBegin(diagonal));
	PW(VecAssemblyEnd(diagonal));

	for (PetscInt q = 1; q < num_qbits + 1; q++) {
		PetscInt side_neighbor = q + 1;
		if (q % num_diagonal_qbits == 0) {
			side_neighbor -= num_diagonal_qbits;
		}
		PetscInt down_neighbor = q + num_diagonal_qbits;
		if (down_neighbor > num_qbits) {
			down_neighbor = down_neighbor % num_qbits;
		}

		if (side_neighbor != q - 1)
			generate_adjacent_interaction_petsc(q, side_neighbor, &diagonal);
		if (down_neighbor != q - num_diagonal_qbits)
			generate_adjacent_interaction_petsc(q, down_neighbor, &diagonal);
	}

	int vecsize, matsize;
	VecGetSize(diagonal, &vecsize);

	PW(MatDiagonalSet(J, diagonal, INSERT_VALUES));
	PW(VecDestroy(&diagonal));

	PW(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
	PW(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
}

void IsingHamiltonian::generate_j_1d() {
	PW(MatCreate(PETSC_COMM_WORLD, &J));
	PW(MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, num_diagonal_terms, num_diagonal_terms));
	PW(MatSetFromOptions(J));
	PW(MatMPIAIJSetPreallocation(J, 1, NULL, 0, NULL));

	Vec diagonal;
	PW(VecCreate(PETSC_COMM_WORLD, &diagonal));
	PW(VecSetSizes(diagonal, PETSC_DECIDE, num_diagonal_terms));
	PW(VecSetFromOptions(diagonal));
	PW(VecZeroEntries(diagonal));
	PW(VecAssemblyBegin(diagonal));
	PW(VecAssemblyEnd(diagonal));

	for (PetscInt q = 1; q < num_qbits + 1; q++) {
		PetscInt side_neighbor = q + 1;
		if (q == num_qbits) {
			side_neighbor = 1;
		}
		generate_adjacent_interaction_petsc(q, side_neighbor, &diagonal);
	}

	int vecsize, matsize;
	VecGetSize(diagonal, &vecsize);

	PW(MatDiagonalSet(J, diagonal, INSERT_VALUES));
	PW(VecDestroy(&diagonal));

	PW(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
	PW(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));


}

void IsingHamiltonian::generate_bz() {

	PW(MatCreate(PETSC_COMM_WORLD, &Bz));
	PW(MatSetSizes(Bz, PETSC_DECIDE, PETSC_DECIDE, num_diagonal_terms, num_diagonal_terms));
	PW(MatSetFromOptions(Bz));
	PW(MatSetUp(Bz));

	Vec v;
	PW(VecCreate(PETSC_COMM_WORLD, &v));
	PW(VecSetSizes(v, PETSC_DECIDE, num_diagonal_terms));
	PW(VecSetFromOptions(v));
	PW(VecAssemblyBegin(v));
	PW(VecAssemblyEnd(v));

	MPI_Barrier(MPI_COMM_WORLD);
	PetscInt     local_processor_vector_begin;
	PetscInt     local_processor_vector_end;
	PetscScalar *local_processor_array_pointer;
	PW(VecGetOwnershipRange(v, &local_processor_vector_begin, &local_processor_vector_end));
	PW(VecGetArray(v, &local_processor_array_pointer));

	for (PetscInt x = local_processor_vector_begin; x < local_processor_vector_end; x++) {

		PetscScalar output = 0.0;
		// which qubits are going to be -1 at this location?
		for (PetscInt y = 0; y < num_qbits; y++) {
			output += ((x % (1 << (y + 1))) < (1 << y)) ? 1 : -1;
		}
		local_processor_array_pointer[x - local_processor_vector_begin] = output;
	}

	PW(VecRestoreArray(v, &local_processor_array_pointer));
	MPI_Barrier(MPI_COMM_WORLD);
	PW(MatDiagonalSet(Bz, v, INSERT_VALUES));
	PW(VecDestroy(&v));
	PW(MatAssemblyBegin(Bz, MAT_FINAL_ASSEMBLY));
	PW(MatAssemblyEnd(Bz, MAT_FINAL_ASSEMBLY));
}

void IsingHamiltonian::generate_bx() {

	PW(MatCreate(PETSC_COMM_WORLD, &Bx));
	PW(MatSetSizes(Bx, PETSC_DECIDE, PETSC_DECIDE, num_diagonal_terms, num_diagonal_terms));
	PW(MatSetFromOptions(Bx));
	PW(MatMPIAIJSetPreallocation(Bx, num_qbits, NULL, num_qbits, NULL));

	// Precompute the pattern of the matrix
	int              x = 1;
	std::vector<int> powers_of_two;
	while (x <= num_diagonal_terms) {
		powers_of_two.push_back(x);
		x *= 2;
	}
	// Compute each row from that pattern
	PetscInt first_row_on_processor;
	PetscInt exclusive_row_bound;

	PW(MatGetOwnershipRange(Bx, &first_row_on_processor, &exclusive_row_bound));
	MPI_Barrier(MPI_COMM_WORLD);

	for (PetscInt current_row = first_row_on_processor; current_row < exclusive_row_bound; current_row++) {
		std::vector<PetscInt> colIndex;
		std::vector<double>   row;

		// The element at (row, row) is a reflection point.. that is, the pattern reflects over this point, although this point should be zero
		for (int pattern_idx = 0; pattern_idx < powers_of_two.size(); pattern_idx++) {
			int current_power_of_two = powers_of_two[pattern_idx];
			int column_1             = current_row - current_power_of_two;
			int column_2             = current_row + current_power_of_two;
			if (column_1 >= 0) {
				PetscInt value = ((((column_1 / (current_row - column_1)) + 1) % 2) == 0) ? 0 : 1;
				if (value == 1) {
					row.push_back(1);
					colIndex.push_back(column_1);
				}
			}
			if (column_2 < num_diagonal_terms) {
				PetscInt value = ((((current_row / (column_2 - current_row)) + 1) % 2) == 0) ? 0 : 1;
				if (value == 1) {
					row.push_back(1);
					colIndex.push_back(column_2);
				}
			}
		}
		PetscInt rows_to_set[] = {current_row};
		MatSetValues(Bx, 1, &rows_to_set[0], colIndex.size(), colIndex.data(), row.data(), INSERT_VALUES);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	PW(MatAssemblyBegin(Bx, MAT_FINAL_ASSEMBLY));
	PW(MatAssemblyEnd(Bx, MAT_FINAL_ASSEMBLY));
}

void IsingHamiltonian::save_matrix(Mat *matrix, std::string filename, file_save_mode mode) {
	PetscViewer saver;

	if (mode == PETSC_ASCII) {
		PW(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &saver));
	} else {
		PW(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &saver));
	}
	MatView(*matrix, saver);
	PW(PetscViewerDestroy(&saver));
}

void IsingHamiltonian::write_local_parts_to_file(std::string filename, PetscScalar *vector_pointer, PetscInt local_vector_begin, PetscInt local_vector_end) {
	std::ofstream outputFile(filename, std::ios::app | std::ios::binary);
	for (int i = local_vector_begin; i < local_vector_end; i++) {
		outputFile.write((char *)&vector_pointer[i - local_vector_begin], sizeof(double));
	}
	outputFile.flush();
	outputFile.close();
}

void IsingHamiltonian::sync_save_eigenvector(std::string filename, Vec *vector, double &eigenvalue, double j_val, double bx_val, double bz_val) {
	// Get MPI Information
	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	PetscScalar *vector_pointer;
	PetscInt     local_vector_begin, local_vector_end;
	PW(VecGetOwnershipRange(*vector, &local_vector_begin, &local_vector_end));
	PW(VecGetArray(*vector, &vector_pointer));

	if (rank == 0) {
		std::ofstream outputFile(filename, std::ios::app | std::ios::binary);
		int           size = (num_diagonal_terms * 8) + (4 * 8) + 4;
		outputFile.write((char *)&size, sizeof(int));
		outputFile.write((char *)&eigenvalue, sizeof(double));
		outputFile.write((char *)&j_val, sizeof(double));
		outputFile.write((char *)&bx_val, sizeof(double));
		outputFile.write((char *)&bz_val, sizeof(double));
		outputFile.close();
		write_local_parts_to_file(filename, vector_pointer, local_vector_begin, local_vector_end);
	} else {
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		int buff = 1;
		for (int i = 1; i < num_proc; i++) {
			// Allow rank i to start saving.
			MPI_Recv(&buff, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// Look for acknowlegement when finished saving.
			MPI_Send(&buff, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	} else {
		int buff = 0;
		MPI_Send(&buff, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		write_local_parts_to_file(filename, vector_pointer, local_vector_begin, local_vector_end);
		MPI_Recv(&buff, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	PW(VecRestoreArray(*vector, &vector_pointer));
}
