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
Loss<T>::Loss(std::string loss, 
		T lagrange_multiplier, 
		T trade_off_parameter,
		T random_domain_bound)
	: 
		m_lagrange_multiplier(lagrange_multiplier), 
		m_trade_off_parameter(trade_off_parameter), 
		m_random_domain_bound(random_domain_bound)
{
	set_compute_pointer(loss);
}


template <typename T>
void Loss<T>::quadratic(Network<T>& n, const Batch<T>& target_batch)
{
	n.layers.back().errors = n.layers.back().states - target_batch;
	n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();
}

//template <typename T>
//void Loss<T>::quadratic_plus_schrodinger(Network<T>& n, const Batch<T>& target_batch)
//{
	//int batch_shifted_instance;
	//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		//c2_error(d.target_length(), d.num_training_batches());

	//Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(d.target_length());
	//int eigen_index;
	
	//T E;

	//// calculate physics based cost from derivative w.r.t PSI
	//for (int instance = 0; instance < d.num_training_batches(); ++instance)
	//{
		//batch_shifted_instance = instance + batch * d.num_training_batches();
		//eigen_index = static_cast<int>(instance);

		//// tmp = H*PSI
		//tmp = d.sparse_ham_times_vec(batch_shifted_instance, 
				//n.layers.back().states.col(eigen_index));

		//E = d.energy(batch_shifted_instance);

		//c2_error.col(eigen_index) = 
			//d.sparse_ham_times_vec(batch_shifted_instance, tmp) +
			//n.layers.back().states.col(eigen_index) * E * E - 2 * tmp * E;
	//}
	
	//// grab quadratic errors 
	//n.layers.back().errors = n.layers.back().states - d.training_target_batch(batch);

	//// add to quadratic errors the matrix derivative part of the schrodinger equation
	//n.layers.back().errors += m_lagrange_multiplier * c2_error; 

	//// add chain rule effects
	//n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

//}

//template <typename T>
//void Loss<T>::physics_perturbed_quadratic(Network<T>& n, const Dataset<T>& d, int batch)
//{
	//int batch_shifted_instance;
	//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		//perturbation(d.target_length(), d.num_training_batches());

	//int eigen_index;

	//// calculate schrodinger error for each ground state in batch
	//for (int instance = 0; instance < d.num_training_batches(); ++instance)
	//{
		//batch_shifted_instance = instance + batch * d.num_training_batches();
		//eigen_index = static_cast<int>(instance);
		
		//perturbation.col(eigen_index) =
			//d.sparse_ham_times_vec(batch_shifted_instance, n.layers.back().states.col(eigen_index))-
			//d.energy(batch_shifted_instance) * n.layers.back().states.col(eigen_index);
	//}
	
	//// grab quadratic errors 
	//n.layers.back().errors = n.layers.back().states - d.training_target_batch(batch);

	//// add chain rule effects
	//n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

	//// add perturbation outside normal gradient term
	//n.layers.back().errors += m_trade_off_parameter * perturbation;
//}

//template <typename T>
//void Loss<T>::randomly_perturbed_quadratic(Network<T>& n, const Dataset<T>& d, int batch)
//{
	//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> perturbation(d.target_length(), d.num_training_batches());

	//perturbation.setRandom();
	
	//// grab quadratic errors 
	//n.layers.back().errors = n.layers.back().states - d.training_target_batch(batch);

	//// add chain rule effects
	//n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

	//// add perturbation
	//n.layers.back().errors += m_random_domain_bound * perturbation;
//}

template <typename T>
void Loss<T>::set_lagrange_multiplier(T value)
{
	m_lagrange_multiplier = value;
}

template <typename T>
void Loss<T>::set_trade_off_parameter(T value)
{
	m_trade_off_parameter = value;
}

template <typename T>
void Loss<T>::set_random_domain_bound(T value)
{
	m_random_domain_bound = value;
}

template <typename T>
void Loss<T>::set_compute_pointer(std::string loss)
{
	using namespace std::placeholders;
	if (loss == "bb")
	{
		compute = std::bind(&Loss<T>::quadratic, this, _1, _2);
	}
	//else if (loss == "c2")
	//{
		//compute = std::bind(&Loss<T>::quadratic_plus_schrodinger, this, _1, _2);
	//}
	//else if (loss == "pg")
	//{
		//compute = std::bind(&Loss<T>::physics_perturbed_quadratic, this, _1, _2);
	//}
	//else if (loss == "rd")
	//{
		//compute = std::bind(&Loss<T>::randomly_perturbed_quadratic, this, _1, _2);
	//}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: UNSUPPORTED LOSS FUNCTION\n\n";
		close_ascii_escape();	
		exit(-1);
	}
}

template class Loss<float>;
template class Loss<double>;
