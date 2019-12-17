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

template<typename T>
using Pauli = Eigen::Matrix<T, 2, 2>;

template <typename T>
class Operators  
{
	public:
		Operators(int num_qubits);

		Eigen::SparseMatrix<T> hamiltonian(const Fields<T>& fields) const;
		Eigen::SparseMatrix<T> magnetization(void) const;

	private:
		int num_qubits;
		Pauli<T> sigmax, sigmaz, eye2;

		Eigen::SparseMatrix<T> longitudinal(int num_qubits, int site) const;
		Eigen::SparseMatrix<T> transverse(int num_qubits, int site) const;
		Eigen::SparseMatrix<T> coupling(int num_qubits, int site) const;

		std::vector<Eigen::SparseMatrix<T>> coupling_stencils;
		std::vector<Eigen::SparseMatrix<T>> transverse_stencils;
		std::vector<Eigen::SparseMatrix<T>> longitudinal_stencils;
};


	
#endif /* Operators_hpp */

