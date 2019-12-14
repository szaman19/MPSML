///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Operators.cxx ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Operators.hpp"

Operators::Operators(int num_qubits) : num_qubits(num_qubits)
{
	sigmax << 0, 1, 1, 0;
	sigmaz << 1, 0, 0, -1;
	eye2 = Eigen::Matrix2d::Identity();

	coupling_stencil = coupling(num_qubits, 0);
	transverse_stencil = transverse(num_qubits, 0);
	longitudinal_stencil = longitudinal(num_qubits, 0);

	for (int site = 1; site < num_qubits; ++site)
	{
		coupling_stencil += coupling(num_qubits, site);
		transverse_stencil += transverse(num_qubits, site);
		longitudinal_stencil += longitudinal(num_qubits, site);
	}	
}

Eigen::MatrixXd Operators::hamiltonian(const Fields& fields) const
{
	Eigen::MatrixXd term4, out;

	out = 
		-(fields.coupling / num_qubits) * coupling_stencil + 
		-(fields.transverse / num_qubits) * transverse_stencil + 
		-(fields.longitudinal / num_qubits) * longitudinal_stencil;

	// weak logic for doubles, but errs on the side or correctness
	if (fields.disorder_strength != 0)
	{
		term4 = fields.local_fields[0] * longitudinal(num_qubits, 0);
		for (int site = 1; site < num_qubits; ++site)
			fields.local_fields[site] * longitudinal(num_qubits, site);

		out += term4;	
	}

	return out;
}

Eigen::MatrixXd Operators::longitudinal(int num_qubits, int site) const
{
	std::vector<Eigen::Matrix2d> ops(num_qubits);

	for (auto& i : ops) i = eye2;

	ops[site] = sigmaz;

	Eigen::MatrixXd tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	return tmp;
}

Eigen::MatrixXd Operators::transverse(int num_qubits, int site) const
{
	std::vector<Eigen::Matrix2d> ops(num_qubits);

	for (auto& i : ops) i = eye2;

	ops[site] = sigmax;

	Eigen::MatrixXd tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	return tmp;
}

Eigen::MatrixXd Operators::coupling(int num_qubits, int site) const
{
	if (site == num_qubits - 1)
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, 0);
	else 
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, site+1);
}


Eigen::MatrixXd Operators::magnetization(void) const
{
	return longitudinal_stencil;
}











