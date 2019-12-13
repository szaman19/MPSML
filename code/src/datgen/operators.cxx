///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** operators.cxx ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "operators.hpp"

Eigen::MatrixXd longitudinal(int num_qubits, int site)
{
	std::vector<Eigen::Matrix2d> ops(num_qubits);

	for (auto& i : ops) i = eye();

	ops[site] = sigmaz();

	Eigen::MatrixXd tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	return tmp;
}

Eigen::MatrixXd transverse(int num_qubits, int site)
{
	std::vector<Eigen::Matrix2d> ops(num_qubits);

	for (auto& i : ops) i = eye();

	ops[site] = sigmax();

	Eigen::MatrixXd tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	return tmp;
}

Eigen::MatrixXd coupling(int num_qubits, int site)
{
	if (site == num_qubits - 1)
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, 0);
	else 
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, site+1);
}

Eigen::MatrixXd hamiltonian(int num_qubits, const Fields& fields)
{
	Eigen::MatrixXd term1 = coupling(num_qubits, 0);
	Eigen::MatrixXd term2 = transverse(num_qubits, 0);
	Eigen::MatrixXd term3 = longitudinal(num_qubits, 0);
	Eigen::MatrixXd term4 = fields.local_fields[0] * longitudinal(num_qubits, 0);

	for (int site = 1; site < num_qubits; ++site)
	{
		term1 += coupling(num_qubits, site);
		term2 += transverse(num_qubits, site);
		term3 += longitudinal(num_qubits, site);
		term4 += fields.local_fields[site] * longitudinal(num_qubits, site);
	}	

	return 
		-(fields.coupling / num_qubits) * term1 + 
		-(fields.transverse / num_qubits) * term2 +
		-(fields.longitudinal / num_qubits) * term3 + term4;
}

Eigen::VectorXd magnetization(int num_qubits)
{
	Eigen::MatrixXd tmp = longitudinal(num_qubits, 0);

	for (int site = 1; site < num_qubits; ++site)
		tmp += longitudinal(num_qubits, site);

	return tmp.diagonal();
}











