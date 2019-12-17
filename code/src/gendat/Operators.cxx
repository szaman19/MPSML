///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Operators.cxx ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <gendat/Operators.hpp>

template <typename T>
Operators<T>::Operators(int num_qubits) : num_qubits(num_qubits)
{
	sigmax << 0, 1, 1, 0;
	sigmaz << 1, 0, 0, -1;
	eye2 = Pauli<T>::Identity();

	for (int site = 0; site < num_qubits; ++site)
	{
		coupling_stencils.push_back(coupling(num_qubits, site));
		transverse_stencils.push_back(transverse(num_qubits, site));
		longitudinal_stencils.push_back(longitudinal(num_qubits, site));
	}	

	for (int site = 0; site < num_qubits; ++site)
	{
		coupling_stencils[site].makeCompressed();
		longitudinal_stencils[site].makeCompressed();
		transverse_stencils[site].makeCompressed();
	}
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::hamiltonian(const Fields<T>& fields) const
{
	Eigen::SparseMatrix<T> term1, term2, term3;

	term1 = fields.coupling[0] * coupling_stencils[0];
	term2 = fields.transverse[0] * transverse_stencils[0];
	term3 = fields.longitudinal[0] * longitudinal_stencils[0];

	for (int qubit = 1; qubit < num_qubits; ++qubit)
	{
		term1 += fields.coupling[qubit] * coupling_stencils[qubit];
		term2 += fields.transverse[qubit] * transverse_stencils[qubit];
		term3 += fields.longitudinal[qubit] * longitudinal_stencils[qubit];
	}

	return (-term1 - term2 - term3) / num_qubits;
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::longitudinal(int num_qubits, int site) const
{
	std::vector<Pauli<T>> ops(num_qubits);

	for (auto& i : ops) i = eye2;

	ops[site] = sigmaz;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	Eigen::SparseMatrix<T> out(tmp.rows(), tmp.cols());

	for (int i = 0; i < tmp.cols(); ++i)
	{
		for (int j = 0; j < tmp.rows(); ++j)
		{
			if (tmp(j,i) != 0) out.insert(i,j) = tmp(j,i);
		}
	}
	
	return out;
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::transverse(int num_qubits, int site) const
{
	std::vector<Pauli<T>> ops(num_qubits);

	for (auto& i : ops) i = eye2;

	ops[site] = sigmax;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	Eigen::SparseMatrix<T> out(tmp.rows(), tmp.cols());

	for (int i = 0; i < tmp.cols(); ++i)
	{
		for (int j = 0; j < tmp.rows(); ++j)
		{
			if (tmp(j,i) != 0) out.insert(i,j) = tmp(j,i);
		}
	}

	return out;
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::coupling(int num_qubits, int site) const
{
	if (site == num_qubits - 1)
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, 0);
	else 
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, site+1);
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::magnetization(void) const
{
	Eigen::SparseMatrix<T> tmp = longitudinal_stencils[0];

	for (int i = 1; i < longitudinal_stencils.size(); ++i)
		tmp += longitudinal_stencils[i];

	return tmp;
}

template class Operators<float>;
template class Operators<double>;
