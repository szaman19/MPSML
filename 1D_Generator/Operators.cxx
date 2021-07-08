///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Operators.cxx ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Operators.hpp>

template <typename T>
Operators<T>::Operators(int num_qubits) : num_qubits(num_qubits)
{
	std::complex<T> i(0,1);

	sigmax << 0, 1, 1, 0;
	sigmaz << 1, 0, 0, -1;
	sigmay << 0, i, -i, 0;
	eye2 = Pauli<T>::Identity();

	for (int site = 0; site < num_qubits; ++site)
	{
		x_coupling_stencils.push_back(x_coupling(num_qubits, site));
		y_coupling_stencils.push_back(y_coupling(num_qubits, site));
		z_coupling_stencils.push_back(z_coupling(num_qubits, site));

		transverse_stencils.push_back(transverse(num_qubits, site));
		longitudinal_stencils.push_back(longitudinal(num_qubits, site));
	}	

	for (int site = 0; site < num_qubits; ++site)
	{
		x_coupling_stencils[site].makeCompressed();
		y_coupling_stencils[site].makeCompressed();
		z_coupling_stencils[site].makeCompressed();

		longitudinal_stencils[site].makeCompressed();
		transverse_stencils[site].makeCompressed();
	}
}

template <typename T>
int Operators<T>::ham_nnz(std::string model)
{
	Eigen::SparseMatrix<T> term1, term2, term3, out;
		
	term1 = z_coupling_stencils[0];
	term2 = transverse_stencils[0];
	term3 = longitudinal_stencils[0];

	for (int qubit = 1; qubit < num_qubits; ++qubit)
	{
		term1 += z_coupling_stencils[qubit];
		term2 += transverse_stencils[qubit];
		term3 += longitudinal_stencils[qubit];
	}

	out = term1 + term2 + term3;

	return out.nonZeros();

}

template <typename T> Eigen::SparseMatrix<T> 
Operators<T>::ising_hamiltonian(const T* ptr_to_instance_of_fields_batch) const
{
	Eigen::SparseMatrix<T> term1, term2, term3;

	const T* ptr = ptr_to_instance_of_fields_batch;

	term1 = *(ptr + 0*num_qubits) * z_coupling_stencils[0];
	term2 = *(ptr + 1*num_qubits) * transverse_stencils[0];
	term3 = *(ptr + 2*num_qubits) * longitudinal_stencils[0];

	for (int qubit = 1; qubit < num_qubits; ++qubit)
	{
		term1 += *(ptr + qubit + 0*num_qubits) * z_coupling_stencils[qubit];
		term2 += *(ptr + qubit + 1*num_qubits) * transverse_stencils[qubit];
		term3 += *(ptr + qubit + 2*num_qubits) * longitudinal_stencils[qubit];
	}
        
        
	return (-term1 - term2 - term3) / num_qubits;
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::ising_hamiltonian(const Fields<T>& fields) const
{
	Eigen::SparseMatrix<T> term1, term2, term3;

	term1 = fields.coupling[0] * z_coupling_stencils[0];
	term2 = fields.transverse[0] * transverse_stencils[0];
	term3 = fields.longitudinal[0] * longitudinal_stencils[0];

	for (int qubit = 1; qubit < num_qubits; ++qubit)
	{
		term1 += fields.coupling[qubit] * z_coupling_stencils[qubit];
		term2 += fields.transverse[qubit] * transverse_stencils[qubit];
		term3 += fields.longitudinal[qubit] * longitudinal_stencils[qubit];
	}
        
        //std::cout << Eigen::MatrixXd((-term1-term2-term3)/num_qubits) << std::endl;
	return (-term1 - term2 - term3) / num_qubits;
}

// template <typename T> Eigen::SparseMatrix<T> 
// Operators<T>::xyz_hamiltonian(const T* ptr_to_instance_of_fields_batch) const
// {
// 	Eigen::SparseMatrix<T> term1, term2, term3;

// 	const T* ptr = ptr_to_instance_of_fields_batch;

// 	term1 = (XYZ_J + *(ptr + 0*num_qubits)) * x_coupling_stencils[0] + 
// 		    (XYZ_J - *(ptr + 0*num_qubits)) * y_coupling_stencils[0];

// 	term2 = *(ptr + 1*num_qubits) * z_coupling_stencils[0];
// 	term3 = *(ptr + 2*num_qubits) * longitudinal_stencils[0];

// 	for (int qubit = 1; qubit < num_qubits; ++qubit)
// 	{
// 		term1 += (XYZ_J + *(ptr + qubit + 0*num_qubits)) * x_coupling_stencils[qubit] + 
// 			     (XYZ_J - *(ptr + qubit + 0*num_qubits)) * y_coupling_stencils[qubit];

// 		term2 += *(ptr + qubit + 1*num_qubits) * z_coupling_stencils[qubit];

// 		term3 += *(ptr + qubit + 2*num_qubits) * longitudinal_stencils[qubit] * 
// 			std::pow(-1, qubit);
// 	}

// 	return -0.5 * (term1 + term2 + term3) / num_qubits;
// }

// template <typename T> Eigen::SparseMatrix<T> 
// Operators<T>::xyz_hamiltonian(const Fields<T>& fields) const
// {
// 	Eigen::SparseMatrix<T> term1, term2, term3;

// 	T lambda, delta, h;

// 	lambda = fields.coupling[0];
// 	delta = fields.transverse[0]; 
// 	h = fields.longitudinal[0];

// 	term1 = (XYZ_J + lambda) * x_coupling_stencils[0] + 
// 			(XYZ_J - lambda) * y_coupling_stencils[0];

// 	term2 = delta * z_coupling_stencils[0];
// 	term3 = h * longitudinal_stencils[0];

// 	for (int qubit = 1; qubit < num_qubits; ++qubit)
// 	{
// 		lambda = fields.coupling[qubit];
// 		delta = fields.transverse[qubit]; 
// 		h = fields.longitudinal[qubit];

// 		term1 = (XYZ_J + lambda) * x_coupling_stencils[qubit] + 
// 				(XYZ_J - lambda) * y_coupling_stencils[qubit];

// 		term2 = delta * z_coupling_stencils[qubit];
// 		term3 = h * longitudinal_stencils[qubit] * std::pow(-1, qubit);
// 	}

// 	return -0.5 * (term1 + term2 + term3) / num_qubits;
// }
template <typename T>
Eigen::SparseMatrix<T> Operators<T>::energy(const T* ptr_to_instance_of_energy_batch) const
{
	int dim = z_coupling_stencils[0].rows();
	Eigen::SparseMatrix<T> out(dim, dim);

	for (int i = 0; i < dim; ++i)
		out.insert(i, i) = *(ptr_to_instance_of_energy_batch + i);

	out.makeCompressed();

	return out;
}

template <typename T>
Eigen::SparseMatrix<std::complex<T>> Operators<T>::surface(int num_qubits, int site) const
{
	std::vector< Eigen::Matrix<std::complex<T>, 2, 2> > ops(num_qubits);

	for (auto& i : ops) i = eye2;

	ops[site] = sigmay;

	Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> tmp = ops.back();

	for (int i = ops.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(ops[i], tmp).eval();

	Eigen::SparseMatrix<std::complex<T>> out(tmp.rows(), tmp.cols());

	for (int i = 0; i < tmp.cols(); ++i)
	{
		for (int j = 0; j < tmp.rows(); ++j)
		{
			if (tmp(j,i).real() != 0 || tmp(j,i).imag() != 0) out.insert(i,j) = tmp(j,i);
		}
	}
	
	return out;
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
Eigen::SparseMatrix<T> Operators<T>::x_coupling(int num_qubits, int site) const
{
	if (site == num_qubits - 1)
		return transverse(num_qubits, site) * transverse(num_qubits, 0);
	else 
		return transverse(num_qubits, site) * transverse(num_qubits, site+1);
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::y_coupling(int num_qubits, int site) const
{
	if (site == num_qubits - 1)
		return (surface(num_qubits, site) * surface(num_qubits, 0)).real();
	else 
		return (surface(num_qubits, site) * surface(num_qubits, site+1)).real();
}

template <typename T>
Eigen::SparseMatrix<T> Operators<T>::z_coupling(int num_qubits, int site) const
{
	if (site == num_qubits - 1)
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, 0);
	else 
		return longitudinal(num_qubits, site) * longitudinal(num_qubits, site+1);
}


template class Operators<double>;
template class Operators<float>;
