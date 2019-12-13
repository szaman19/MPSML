///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** operators.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef operators_hpp
#define operators_hpp

#include <Eigen/KroneckerProduct>
#include "pauli.hpp"
#include "Fields.hpp"


Eigen::MatrixXd longitudinal(int num_qubits, int site);

Eigen::MatrixXd transverse(int num_qubits, int site);

Eigen::MatrixXd coupling(int num_qubits, int site);

Eigen::MatrixXd hamiltonian(int num_qubits, const Fields& fields);

Eigen::VectorXd magnetization(int num_qubits);
	
#endif /* operators_hpp */

