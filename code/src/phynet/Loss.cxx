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
void Loss<T>::quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
	{
		nets[vec].layers.back().errors = 
			nets[vec].layers.back().states - data.training_target_batch(batch, vec);

		nets[vec].layers.back().errors.array() *= 
			nets[vec].layers.back().derivative_of_activation_on_weighted_sum();
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
	for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
	{
		nets[vec].layers.back().errors = 
			nets[vec].layers.back().states - data.training_target_batch(batch, vec);
	}

	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), SCHRO(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		// 2. 
		for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
			P.col(vec) = nets[vec].layers.back().states.col(instance);

		// 3. 
		H = this->operators.hamiltonian(data.training_feature_batch(batch).col(instance).data());

		// 4. 
		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		// 5. Schrodinger loss
		SCHRO = H * H * P * L + P * E * L * E - (H * P * E * L + H * P * L * E);

		// 6. 
		for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
			nets[vec].layers.back().errors.col(instance) += (1.0 / dim) * SCHRO.col(instance);
	}	

	// 7.
	for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
	{
		nets[vec].layers.back().errors.array() *= 
			nets[vec].layers.back().derivative_of_activation_on_weighted_sum();
	}
}

template <typename T>
void Loss<T>::physics_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
	{
		nets[vec].layers.back().errors = 
			nets[vec].layers.back().states - data.training_target_batch(batch, vec);
		
		nets[vec].layers.back().errors.array() *= 
			nets[vec].layers.back().derivative_of_activation_on_weighted_sum();
	}
	
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P(dim, dim), PERT(dim, dim);
	Eigen::SparseMatrix<T> H(dim, dim), E(dim, dim), L(lagrange_matrix);

	for (int instance = 0; instance < data.batch_size; ++instance)
	{
		for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
			P.col(vec) = nets[vec].layers.back().states.col(instance);

		H = this->operators.hamiltonian(data.training_feature_batch(batch).col(instance).data());

		E = this->operators.energy(data.training_energy_batch(batch).col(instance).data());	

		// 5. Schrodinger perturbation 
		PERT = trade_off_parameter * (H * P - P * E);

		for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
			nets[vec].layers.back().errors.col(instance) += (1.0 / dim) * PERT.col(instance);
	}	
}

template <typename T>
void Loss<T>::randomly_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch)
{
	int dim = data.num_eigenvectors(); 

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> perturbation(dim, data.batch_size);

	perturbation.setRandom();
	
	for (int vec = 0; vec < data.num_eigenvectors(); ++vec)
	{
		nets[vec].layers.back().errors = 
			nets[vec].layers.back().states - data.training_target_batch(batch, vec);
		
		nets[vec].layers.back().errors.array() *= 
			nets[vec].layers.back().derivative_of_activation_on_weighted_sum();

		nets[vec].layers.back().errors += random_domain_bound * perturbation;
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
	else if (loss == "pg")
	{
		compute = std::bind(&Loss<T>::physics_perturbed_quadratic, this, _1, _2, _3);
	}
	else if (loss == "rd")
	{
		compute = std::bind(&Loss<T>::randomly_perturbed_quadratic, this, _1, _2, _3);
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
