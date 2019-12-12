///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Model.hpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Model.hpp"
		
Model::Model(const Network& net, const Loss& loss, const Optimizer& optimizer) 
	: m_net(net), m_loss(loss), m_optimizer(optimizer) 
{ 
 	// Null constructor body 
}

void Model::learn_from(const Dataset& dataset)
{
	for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
	{
		m_net.feedforward(dataset.feature_batch(batch));
		m_loss.compute(m_net, dataset, batch);
		m_net.backpropagate();
		m_optimizer.update(m_net);
	}
}

double Model::mse(const Dataset& dataset)
{
	Eigen::ArrayXXd tmp(dataset.target_length(), dataset.batch_size());
	tmp.setZero();

	for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
	{
		m_net.feedforward(dataset.feature_batch(batch));
		m_loss.compute(m_net, dataset, batch);
		tmp += (m_net.layers.back().errors).array().square();
	}

	return tmp.mean() / dataset.batches();

}

void Model::write_schrodinger_error(const Dataset &dataset, std::string filename, int epoch)
{
	std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");
	long idx;

	std::size_t batch_shifted_instance;
	Eigen::MatrixXd c2_error(dataset.target_length(), dataset.batch_size());
	Eigen::VectorXd tmp(dataset.target_length());
	int eigen_index;
	double E;

	if (file.is_open())
	{
		file << "#Field   ||Schrodinger Error||^2   <Schrodinger Error>  \n";
		file << std::scientific;

		for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
		{
			m_net.feedforward(dataset.feature_batch(batch));

			for (std::size_t instance = 0; instance < dataset.batch_size(); ++instance)
			{
				idx = static_cast<long>(instance);

				for (std::size_t i = 0; i < dataset.feature_length(); ++i)
				{
					file << std::setw(8);
					file << m_net.layers[0].states.col(idx)(static_cast<long>(i)) << '\t';
				}


				batch_shifted_instance = instance + batch * dataset.batch_size();
				eigen_index = static_cast<std::size_t>(instance);

				tmp = dataset.sparse_ham_times_vec(batch_shifted_instance, 
						m_net.layers.back().states.col(eigen_index));

				E = dataset.energy(batch_shifted_instance);

				c2_error.col(eigen_index) = 
					dataset.sparse_ham_times_vec(batch_shifted_instance, tmp) +
					m_net.layers.back().states.col(eigen_index) * E * E - 2 * tmp * E;

				file << std::setw(8) << std::right;
				file << c2_error.norm() * c2_error.norm() << '\t';
				file << std::setw(8) << std::right;
				file << c2_error.mean() << '\n';
				
			}
		}
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO WRITE SCHRODINGER ERROR\n";
		std::cerr << std::endl;
		close_ascii_escape();
	}

}


void Model::write_average_magnetization(const Dataset &dataset, std::string filename, int epoch)
{
	std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");
	long idx, idi;

	double avg_of_value, avg_of_square;
	double coeff_squared, coeff;

	Eigen::VectorXd szdiag = dataset.szdiag();

	if (file.is_open())
	{
		file << "#Field   <Sz>   <Sz^2>   <Sz>^2 \n";
		file << std::scientific;

		for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
		{
			m_net.feedforward(dataset.feature_batch(batch));

			for (std::size_t instance = 0; instance < dataset.batch_size(); ++instance)
			{
				idx = static_cast<long>(instance);

				for (std::size_t i = 0; i < dataset.feature_length(); ++i)
				{
					file << std::setw(8);
					file << m_net.layers[0].states.col(idx)(static_cast<long>(i)) << '\t';
				}

				avg_of_value = 0;
				avg_of_square = 0;

				for (std::size_t i = 0; i < dataset.target_length(); ++i) 
				{
					idi = static_cast<long>(i);
					coeff = m_net.layers.back().states.col(idx)(idi);
					coeff_squared = coeff * coeff;
					avg_of_value += szdiag(idi) * coeff_squared; 
					avg_of_square += szdiag(idi) * szdiag(idi) * coeff_squared;
				}
				
				file << std::setw(8) << std::right;
				file << avg_of_value << '\t';
				file << std::setw(8) << std::right;
				file << avg_of_square << '\t';	
				file << std::setw(8) << std::right;
				file << avg_of_value * avg_of_value << '\t';
				file << '\n';
			}
		}
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO WRITE COEFFICIENTS\n";
		std::cerr << std::endl;
		close_ascii_escape();
	}

}

void Model::write_coefficients(const Dataset &dataset, std::string filename, int epoch)
{
	std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");
	long idx;

	if (file.is_open())
	{
		file << std::scientific;

		for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
		{
			m_net.feedforward(dataset.feature_batch(batch));

			for (std::size_t instance = 0; instance < dataset.batch_size(); ++instance)
			{
				idx = static_cast<long>(instance);

				for (std::size_t i = 0; i < dataset.feature_length(); ++i)
				{
					file << std::setw(8);
					file << m_net.layers[0].states.col(idx)(static_cast<long>(i)) << '\t';
				}

				for (std::size_t i = 0; i < dataset.target_length(); ++i) 
				{
					file << std::setw(20) << std::right;
					file << m_net.layers.back().states.col(idx)(static_cast<long>(i)) << '\t';
				}
				
				file << '\n';
			}
		}
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO WRITE COEFFICIENTS\n";
		std::cerr << std::endl;
		close_ascii_escape();
	}
}


void Model::write_lyapunov_estimate(const Dataset &dataset, std::string filename, int epoch)
{
	std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");

	Eigen::ArrayXXd perturbation(dataset.feature_length(), dataset.batch_size());
	Eigen::ArrayXXd output_batch(dataset.target_length(), dataset.batch_size());

	double rise, run;
	long idx;
	
	if (file.is_open())
	{
		file << std::scientific;

		for (std::size_t batch = 0; batch < dataset.batches(); ++batch)
		{
			m_net.feedforward(dataset.feature_batch(batch));
			output_batch = m_net.layers.back().states;

			perturbation.setRandom();
			perturbation *= 0.001;

			m_net.feedforward(dataset.feature_batch(batch) + perturbation.matrix());

			for (std::size_t instance = 0; instance < dataset.batch_size(); ++instance)
			{
				idx = static_cast<long>(instance);

				for (std::size_t i = 0; i < dataset.feature_length(); ++i)
				{
					file << std::setw(8);
					file << m_net.layers[0].states.col(idx)(static_cast<long>(i)) << '\t';
				}

				rise = (output_batch.col(idx).matrix() - m_net.layers.back().states.col(idx)).norm();
				run = perturbation.matrix().col(idx).norm();

				file << rise / run << '\n';
			}
		}
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO WRITE LYAPUNOV EXPONENT\n";
		std::cerr << std::endl;
		close_ascii_escape();
	}

}


void Model::save(std::string filename)
{
	std::ofstream file(filename.c_str());

	if (file.is_open())
	{
		std::cout << "Saving net to " << filename << '\n';

		file << "Number of layers:   " << m_net.layers.size() << '\n';

		for (std::size_t i = 0; i < m_net.layers.size(); ++i)
			file << "Neurons in layer " << i << ": " << m_net.layers[i].biases.rows() << '\n';

		file << std::scientific;

		for (std::size_t i = 1; i < m_net.layers.size(); ++i)
		{
			file << "layer: " << i << '\n';
			file << "weights: \n";
			file << m_net.layers[i].weights << '\n';
			file << "biases: \n";
			file << m_net.layers[i].biases.col(0) << '\n';
			file << "--------------\n";
		}	
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO SAVE NET\n";
		std::cerr << std::endl;
		close_ascii_escape();
	}
}

