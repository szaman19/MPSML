///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Operators.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Operators_hpp
#define Operators_hpp

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include "Fields.hpp"

#define XYZ_J 1.0

template<typename T>
using Pauli = Eigen::Matrix<T, 2, 2>;

template <typename T>
class Operators  
{
	public:
		Operators(int num_qubits);

		Eigen::SparseMatrix<T> ising_hamiltonian(const Fields<T>& fields) const;
		Eigen::SparseMatrix<T> ising_hamiltonian(const T* ptr_to_instance_of_fields_batch) const;

		Eigen::SparseMatrix<T> energy(const T* ptr_to_instance_of_energy_batch) const;


		int ham_nnz(std::string model);
		const int num_qubits;

	private:
		Pauli<T> sigmax, sigmaz, eye2;
		Eigen::Matrix<std::complex<T>, 2, 2> sigmay;

		Eigen::SparseMatrix<T> longitudinal(int num_qubits, int site) const;
		Eigen::SparseMatrix<T> transverse(int num_qubits, int site) const;
		Eigen::SparseMatrix<std::complex<T>> surface(int num_qubits, int site) const;

		Eigen::SparseMatrix<T> x_coupling(int num_qubits, int site) const;
		Eigen::SparseMatrix<T> y_coupling(int num_qubits, int site) const;
		Eigen::SparseMatrix<T> z_coupling(int num_qubits, int site) const;

		std::vector<Eigen::SparseMatrix<T>> x_coupling_stencils;
		std::vector<Eigen::SparseMatrix<T>> y_coupling_stencils;
		std::vector<Eigen::SparseMatrix<T>> z_coupling_stencils;

		std::vector<Eigen::SparseMatrix<T>> transverse_stencils;
		std::vector<Eigen::SparseMatrix<T>> surface_stencils;
		std::vector<Eigen::SparseMatrix<T>> longitudinal_stencils;
};


	
#endif /* Operators_hpp */

