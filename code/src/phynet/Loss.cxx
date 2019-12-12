///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Loss.cpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Loss.hpp"

Loss::Loss(std::string loss, 
		double lagrange_multiplier, 
		double trade_off_parameter,
		double random_domain_bound)
	: 
		m_lagrange_multiplier(lagrange_multiplier), 
		m_trade_off_parameter(trade_off_parameter), 
		m_random_domain_bound(random_domain_bound)
{
	set_compute_pointer(loss);
}


void Loss::quadratic(Network& n, const Dataset& d, std::size_t batch)
{
	n.layers.back().errors = n.layers.back().states - d.target_batch(batch);
	n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();
}

void Loss::quadratic_plus_schrodinger(Network& n, const Dataset& d, std::size_t batch)
{
	std::size_t batch_shifted_instance;
	Eigen::MatrixXd c2_error(d.target_length(), d.batch_size());
	Eigen::VectorXd tmp(d.target_length());
	int eigen_index;
	double E;

	// calculate physics based cost from derivative w.r.t PSI
	for (std::size_t instance = 0; instance < d.batch_size(); ++instance)
	{
		batch_shifted_instance = instance + batch * d.batch_size();
		eigen_index = static_cast<std::size_t>(instance);

		// tmp = H*PSI
		tmp = d.sparse_ham_times_vec(batch_shifted_instance, 
				n.layers.back().states.col(eigen_index));

		E = d.energy(batch_shifted_instance);

		c2_error.col(eigen_index) = 
			d.sparse_ham_times_vec(batch_shifted_instance, tmp) +
			n.layers.back().states.col(eigen_index) * E * E - 2 * tmp * E;
	}
	
	// grab quadratic errors 
	n.layers.back().errors = n.layers.back().states - d.target_batch(batch);

	// add to quadratic errors the matrix derivative part of the schrodinger equation
	n.layers.back().errors += m_lagrange_multiplier * c2_error; 

	// add chain rule effects
	n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

}

void Loss::physics_perturbed_quadratic(Network& n, const Dataset& d, std::size_t batch)
{
	std::size_t batch_shifted_instance;
	Eigen::MatrixXd perturbation(d.target_length(), d.batch_size());
	int eigen_index;

	// calculate schrodinger error for each ground state in batch
	for (std::size_t instance = 0; instance < d.batch_size(); ++instance)
	{
		batch_shifted_instance = instance + batch * d.batch_size();
		eigen_index = static_cast<std::size_t>(instance);
		
		perturbation.col(eigen_index) =
			d.sparse_ham_times_vec(batch_shifted_instance, n.layers.back().states.col(eigen_index))-
			d.energy(batch_shifted_instance) * n.layers.back().states.col(eigen_index);
	}
	
	// grab quadratic errors 
	n.layers.back().errors = n.layers.back().states - d.target_batch(batch);

	// add chain rule effects
	n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

	// add perturbation outside normal gradient term
	n.layers.back().errors += m_trade_off_parameter * perturbation;
}

void Loss::randomly_perturbed_quadratic(Network& n, const Dataset& d, std::size_t batch)
{
	Eigen::MatrixXd perturbation(d.target_length(), d.batch_size());

	perturbation.setRandom();
	
	// grab quadratic errors 
	n.layers.back().errors = n.layers.back().states - d.target_batch(batch);

	// add chain rule effects
	n.layers.back().errors.array() *= n.layers.back().derivative_of_activation_on_weighted_sum();

	// add perturbation
	n.layers.back().errors += m_random_domain_bound * perturbation;
}

void Loss::set_lagrange_multiplier(double value)
{
	m_lagrange_multiplier = value;
}

void Loss::set_trade_off_parameter(double value)
{
	m_trade_off_parameter = value;
}

void Loss::set_random_domain_bound(double value)
{
	m_random_domain_bound = value;
}

void Loss::set_compute_pointer(std::string loss)
{
	using namespace std::placeholders;
	if (loss == "bb")
	{
		compute = std::bind(&Loss::quadratic, this, _1, _2, _3);
	}
	else if (loss == "c2")
	{
		compute = std::bind(&Loss::quadratic_plus_schrodinger, this, _1, _2, _3);
	}
	else if (loss == "pg")
	{
		compute = std::bind(&Loss::physics_perturbed_quadratic, this, _1, _2, _3);
	}
	else if (loss == "rd")
	{
		compute = std::bind(&Loss::randomly_perturbed_quadratic, this, _1, _2, _3);
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: UNSUPPORTED LOSS FUNCTION\n\n";
		close_ascii_escape();	
		exit(-1);
	}
}
