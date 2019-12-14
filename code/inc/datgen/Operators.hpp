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

#include <Eigen/KroneckerProduct>
#include "Fields.hpp"

class Operators  
{
	public:
		Operators(int num_qubits);

		Eigen::MatrixXd hamiltonian(const Fields& fields) const;
		Eigen::MatrixXd magnetization(void) const;

	private:
		int num_qubits;
		Eigen::Matrix2d sigmax, sigmaz, eye2;

		Eigen::MatrixXd longitudinal(int num_qubits, int site) const;
		Eigen::MatrixXd transverse(int num_qubits, int site) const;
		Eigen::MatrixXd coupling(int num_qubits, int site) const;

		Eigen::MatrixXd coupling_stencil;
		Eigen::MatrixXd transverse_stencil;
		Eigen::MatrixXd longitudinal_stencil;
};

	
#endif /* Operators_hpp */

