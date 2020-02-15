///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Loss.cpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Loss.hpp>

template <typename T>
Loss<T>::Loss(std::string loss, const Operators<T>& operators,
		T lagrange_multiplier, 
		T trade_off_parameter,
		T random_domain_bound)
	: 
		lagrange_multiplier(lagrange_multiplier), 
		trade_off_parameter(trade_off_parameter), 
		random_domain_bound(random_domain_bound),
		operators(operators)
{
	set_compute_pointer(loss);

	int dim = operators.magnetization().rows();

	lagrange_matrix.resize(dim, dim);

	for (int i = 0; i < dim; ++i)
		lagrange_matrix.insert(i,i) = lagrange_multiplier;

	lagrange_matrix.makeCompressed();
}

template <typename T>
Eigen::RowVectorXd Loss<T>::pure_cost(NetVec<T>& nets, const Dataset<T>& data)
{
	int dim = data.num_eigenvectors(); 
	
	Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor> out(2, dim);

	out.setZero();

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim);

	for (int batch = 0; batch < data.num_testing_batches(); ++batch)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].feedforward(data.testing_feature_batch(batch));

		for (int instance = 0; instance < data.batch_size; ++instance)
		{
			for (std::size_t net = 0; net < nets.size(); ++net)
			{
				P.col(net) = nets[net].layers.back().states.col(instance) / 
					nets[net].layers.back().states.col(instance).norm();
			}

			H = this->operators.ising_hamiltonian(data.testing_feature_batch(batch).col(instance).data());

			E = this->operators.energy(data.testing_energy_batch(batch).col(instance).data());	

			for (std::size_t net = 0; net < nets.size(); ++net)
			{
				out(0, net) += 0.5 * (data.testing_target_batch(batch, net).col(instance) -
					nets[net].layers.back().states.col(instance)).cwiseAbs2().sum();	

				out(1, net) += 0.5 * (H * P - P * E).col(net).cwiseAbs2().sum();
			}
		}	
	}

	Eigen::Map<Eigen::RowVectorXd> v(out.data(), out.size());

	return v / data.num_testing_instances();
}


template <typename T>
void Loss<T>::quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);

		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::quadratic_plus_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	/*
	 1. add the quadratic cost to the errors, can be done in batch mode  
	 2. generate PSI for current instance from states of networks set by Model
	 3. generate sparse H from current instance of inputs 
	 4. generate sparse E for current instance of inputs 
	 5. compute Schrodinger loss as derivative of || (H P - P E) L ||_F.
	 6. take each column of result and add too corresponding column of errors
	 7. add chain rule effect 
	*/

	// 1. 
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
	}

	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		// 2. 
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		// 3. 
		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		// 4. 
		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		// 5. Schrodinger loss
		S = H * H * P * L + P * E * L * E - (H * P * E * L + H * P * L * E);

		// 6. 
		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
		 
	}	

	// 7.
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::sigmoid_semi_supervised(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, 1), S(dim, 1);
	Eigen::SparseMatrix<T> H(dim, dim);

	double E, L, u;

	L = lagrange_matrix.coeff(0,0);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		// if Bx < 0.5, add MSE loss
		if (data.training_feature_batch(batch).col(instance)(data.num_qubits) <= 0.5)
		{
			for (std::size_t net = 0; net < nets.size(); ++net)
			{
				nets[net].layers.back().errors.col(instance) = 
					nets[net].layers.back().states.col(instance) - 
					data.training_target_batch(batch, net).col(instance);
			}
		}

		// 2. 
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		// 3. 
		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		// 4. 
		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data()).coeff(0,0);	

		u = (H * P - E * P).norm();

		// 5. Schrodinger loss
		S = L * sigmoid(0.5*u*u) * (1.0 - sigmoid(0.5*u*u)) * (H * H * P + E * E * P - 2 * E * H * P);

		// 6. 
		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	// 7.
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::supervised_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim), G(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();

			G.col(net) = data.training_target_batch(batch, net).col(instance);
		}

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	
	
		S = (H * H) * (P - G) * L;

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::unsupervised_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		S = H * H * P * L + P * E * L * E - (H * P * E * L + H * P * L * E);

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::shuijun_diannong(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		S = 4 * H * P + 2 * H * P * P.transpose() * P + 2 * P * P.transpose() * H * P;

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::rayleigh_ritz(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::SparseMatrix<T> H(dim, dim);

	Eigen::VectorXd P(dim), S(dim);

	double d, a, l = lagrange_matrix.coeff(0,0);

	if (nets.size() != 1)
	{
		std::cerr << "Only ground state loss\n";
		exit(-10);
	}
	
	nets[0].layers.back().errors = 
		nets[0].layers.back().states - data.training_target_batch(batch, 0);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		P.col(0) = nets[0].layers.back().states.col(instance);

		H = operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		d = P.cwiseAbs2().sum();
		a = P.transpose() * H * P;

		S = 2 * (H * P * d - a * P) / (d*d);

		nets[0].layers.back().errors.col(instance) += l * S;
	}	

	nets[0].layers.back().errors.array() *= 
		nets[0].layers.back().derivative_of_activation_on_weighted_sum();
}

template <typename T>
void Loss<T>::sigmoid_rayleigh_ritz(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::SparseMatrix<T> H(dim, dim);

	Eigen::VectorXd P(dim), S(dim);

	double d, a, l = lagrange_matrix.coeff(0,0);

	if (nets.size() != 1)
	{
		std::cerr << "Only ground state loss\n";
		exit(-10);
	}
	
	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		// if Bx < 0.5, add MSE loss
		if (data.training_feature_batch(batch).col(instance)(data.num_qubits) <= 0.5)
		{
			for (std::size_t net = 0; net < nets.size(); ++net)
			{
				nets[net].layers.back().errors.col(instance) = 
					nets[net].layers.back().states.col(instance) - 
					data.training_target_batch(batch, net).col(instance);
			}
		}

		P.col(0) = nets[0].layers.back().states.col(instance);

		H = operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		d = P.cwiseAbs2().sum();
		a = P.transpose() * H * P;

		S = sigmoid(a/d) * (1 - sigmoid(a/d)) * 2 * (H * P * d - a * P) / (d*d);

		nets[0].layers.back().errors.col(instance) += l * S;
	}	

	nets[0].layers.back().errors.array() *= 
		nets[0].layers.back().derivative_of_activation_on_weighted_sum();
}

template <typename T>
void Loss<T>::raw_sums(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	Eigen::RowVectorXd one(dim), out(dim);

	one.setOnes();

	if (nets.size() != 1)
	{
		std::cerr << "Only ground state loss\n";
		exit(-10);
	}
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
	}

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		P.col(0) = nets[0].layers.back().states.col(instance);

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		out = (E.coeff(0,0) * one * P.col(0) - one * H * P.col(0)) * (E.coeff(0,0) * one - one * H);

		nets[0].layers.back().errors.col(instance) += out.transpose();
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::data_shuijun_diannong(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
	}

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
			P.col(net) = nets[net].layers.back().states.col(instance);

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		S = 4 * H * P + 2 * H * P * P.transpose() * P + 2 * P * P.transpose() * H * P;

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::sigmoid_unitarity(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim), G(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	double u;
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
	}
	
	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance);
				/// nets[net].layers.back().states.col(instance).norm();
			
			G.col(net) = data.training_target_batch(batch, net).col(instance);
		}

		u = std::fabs((P.transpose() * P).trace() - dim);

		S = 4 * sigmoid(u*u) * (1.0 - sigmoid(u*u)) * u * P;

		//std::cout << "u: " << u << '\n';
		//std::cout << S << '\n';
		//exit(-10);

		//if (std::signbit(u) == 1) S *= -1;

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}

}

template <typename T>
void Loss<T>::sigmoid_frobenius(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		P(dim, dim), S(dim, dim), G(dim, dim), I(dim,dim);

	I = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim,dim);

	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	double u;
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
	}
	
	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance);
				/// nets[net].layers.back().states.col(instance).norm();
			
			//G.col(net) = data.training_target_batch(batch, net).col(instance);
		}

		u = (P.transpose() * P - I).norm();

		S = sigmoid(0.25*u*u) * (1.0 - sigmoid(0.25*u*u)) * (P * P.transpose() * P - P);

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}

}

template <typename T>
void Loss<T>::physics_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
		
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}
	
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), PERT(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		// 5. Schrodinger perturbation 
		PERT = trade_off_parameter * (H * P - P * E);

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += PERT.col(instance);
	}	
}

template <typename T>
void Loss<T>::randomly_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> perturbation(dim, data.batch_size);

	perturbation.setRandom();
	
	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors = 
			nets[net].layers.back().states - data.training_target_batch(batch, net);
		
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();

		nets[net].layers.back().errors += random_domain_bound * perturbation;
	}
}

template <typename T>
void Loss<T>::abs_formulation(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), S(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (std::size_t net = 0; net < nets.size(); ++net)
		{
			P.col(net) = nets[net].layers.back().states.col(instance) / 
				nets[net].layers.back().states.col(instance).norm();
		}

		H = this->operators.ising_hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		S = 2 * H * P;

		if (std::signbit(E.coeff(0,0) - (P.transpose() * H * P)(0,0)) == -1) S *= -1;

		for (std::size_t net = 0; net < nets.size(); ++net)
			nets[net].layers.back().errors.col(instance) += S.col(net);
	}	

	for (std::size_t net = 0; net < nets.size(); ++net)
	{
		nets[net].layers.back().errors.array() *= 
			nets[net].layers.back().derivative_of_activation_on_weighted_sum();
	}

}

template <typename T>
void Loss<T>::set_lagrange_multiplier(T value)
{
	lagrange_multiplier = value;
}

template <typename T>
void Loss<T>::set_trade_off_parameter(T value)
{
	trade_off_parameter = value;
}

template <typename T>
void Loss<T>::set_random_domain_bound(T value)
{
	random_domain_bound = value;
}

template <typename T>
void Loss<T>::set_compute_pointer(std::string loss)
{
	using namespace std::placeholders;
	if (loss == "bb")
	{
		compute = std::bind(&Loss<T>::quadratic, this, _1, _2, _3);
	}
	else if (loss == "c2")
	{
		compute = std::bind(&Loss<T>::quadratic_plus_schrodinger, this, _1, _2, _3);
	}
	else if (loss == "un")
	{
		compute = std::bind(&Loss<T>::unsupervised_schrodinger, this, _1, _2, _3);
	}
	else if (loss == "ss")
	{
		compute = std::bind(&Loss<T>::supervised_schrodinger, this, _1, _2, _3);
	}
	else if (loss == "ab")
	{
		compute = std::bind(&Loss<T>::abs_formulation, this, _1, _2, _3);
	}
	else if (loss == "pg")
	{
		compute = std::bind(&Loss<T>::physics_perturbed_quadratic, this, _1, _2, _3);
	}
	else if (loss == "rd")
	{
		compute = std::bind(&Loss<T>::randomly_perturbed_quadratic, this, _1, _2, _3);
	}
	else if (loss == "su")
	{
		compute = std::bind(&Loss<T>::sigmoid_unitarity, this, _1, _2, _3);
	}
	else if (loss == "sd")
	{
		compute = std::bind(&Loss<T>::shuijun_diannong, this, _1, _2, _3);
	}
	else if (loss == "dd")
	{
		compute = std::bind(&Loss<T>::data_shuijun_diannong, this, _1, _2, _3);
	}
	else if (loss == "rs")
	{
		compute = std::bind(&Loss<T>::raw_sums, this, _1, _2, _3);
	}
	else if (loss == "rr")
	{
		compute = std::bind(&Loss<T>::rayleigh_ritz, this, _1, _2, _3);
	}
	else if (loss == "srr")
	{
		compute = std::bind(&Loss<T>::rayleigh_ritz, this, _1, _2, _3);
	}
	else if (loss == "sf")
	{
		compute = std::bind(&Loss<T>::sigmoid_frobenius, this, _1, _2, _3);
	}
	else if (loss == "se")
	{
		compute = std::bind(&Loss<T>::sigmoid_semi_supervised, this, _1, _2, _3);
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: UNSUPPORTED LOSS FUNCTION\n\n";
		close_ascii_escape();	
		exit(-1);
	}
}

template class Loss<double>;
