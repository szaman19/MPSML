///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Model.hpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Model.hpp>

template <typename T>
Model<T>::
Model(const std::vector<Network<T>>& networks, const Loss<T>& loss, const Optimizer<T>& optimizer) 
	: networks(networks), loss(loss), optimizer(optimizer) 
{ 
 	// Null constructor body 
}

template <typename T>
T Model<T>::predictive_power(const Dataset<T>& dataset, int epoch)
{
	int dim = dataset.num_eigenvectors();
	T out = 0;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_preds(dim, dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_targs(dim, dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_final(dim, dim);

	normalized_final.setZero();

	for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
	{
		for (int vec = 0; vec < dataset.num_eigenvectors(); ++vec)
			networks[vec].feedforward(dataset.testing_feature_batch(batch));

		for (int instance = 0; instance < dataset.batch_size; ++instance)
		{
			for (int vec = 0; vec < dataset.num_eigenvectors(); ++vec)
			{
				normalized_targs.col(vec) = 
					dataset.testing_target_batch(batch, vec).col(instance);

				normalized_preds.col(vec) = 
					networks[vec].layers.back().states.col(instance) / 
					networks[vec].layers.back().states.col(instance).norm();
			}

			normalized_final += normalized_preds.transpose() * normalized_targs;

			//std::cout << "#### BATCH: " << batch << " INSTANCE: " << instance << " PREDS";
			//pretty_print(normalized_preds);
			//std::cout << "#### BATCH: " << batch << " INSTANCE: " << instance << " TARGS";
			//pretty_print(normalized_targs);
		}
	}

	normalized_final /= dataset.num_testing_instances();

	std::cout << "Average \n";
	pretty_print(normalized_final);

	return out / (dataset.num_testing_instances() * std::sqrt(dim));	
}

template <typename T>
Eigen::RowVectorXd Model<T>::pure_cost(const Dataset<T>& dataset) 
{
	return loss.pure_cost(networks, dataset);
}

template <typename T>
void Model<T>::pretty_print(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
	std::cout << std::endl;
	std::cout << std::endl;
	for (int i = 0; i < mat.rows(); ++i)
	{
		for (int j = 0; j < mat.cols(); ++j)
		{
			std::cout << std::fixed << std::setprecision(4) << std::setw(9) << std::right;
			if (std::abs(mat(i,j)) == mat.col(j).cwiseAbs().maxCoeff() && i == j) 
			{
				open_ascii_escape("green");
				std::cout << mat(i,j);
				close_ascii_escape();
			}
			else if (std::abs(mat(i,j)) == mat.col(j).cwiseAbs().maxCoeff())
			{
				open_ascii_escape("red");
				std::cout << mat(i,j);
				close_ascii_escape();
			}
			else
			{
				std::cout << mat(i,j);
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;;
}

template <typename T>
void Model<T>::learn_from(const Dataset<T>& dataset)
{
	std::vector<int> batch_indices(dataset.num_training_batches());

	std::iota(batch_indices.begin(), batch_indices.end(), 0);

	std::random_device rd;	
	std::default_random_engine engine(rd());

	std::shuffle(batch_indices.begin(), batch_indices.end(), engine);


	for (auto batch : batch_indices)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (std::size_t net = 0; net < networks.size(); ++net)
			networks[net].feedforward(dataset.training_feature_batch(batch));

		loss.compute(networks, dataset, batch);

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (std::size_t net = 0; net < networks.size(); ++net)
			networks[net].backpropagate();

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (std::size_t net = 0; net < networks.size(); ++net)
			optimizer.update(networks[net]);
	}
}


template <typename T>
T Model<T>::mse(const Dataset<T>& dataset)
{
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> 
		tmp(dataset.target_length(), dataset.batch_size);

	tmp.setZero();

	for (int vec = 0; vec < dataset.num_eigenvectors(); ++vec)
	{
		for (int batch = 0; batch < dataset.num_validation_batches(); ++batch)
		{
			//networks[vec].feedforward(dataset.validation_feature_batch(batch));

			//tmp += (networks[vec].layers.back().states 
					//- dataset.validation_target_batch(batch, vec)).array().square();
			//loss.compute(networks, dataset, batch);

			tmp += networks[vec].layers.back().errors.array().square();
		}
	}
	
	return tmp.mean() / (dataset.num_eigenvectors() * dataset.num_validation_batches());	
}







//template <typename T>
//void Model<T>::write_lyapunov_estimate(const Dataset<T> &dataset, std::string filename, int epoch)
//{
	//std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");

	//Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> 
		//perturbation(dataset.feature_length(), dataset.batch_size);

	//Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> 
		//output_batch(dataset.target_length(), dataset.batch_size);

	//T rise, run;
	//long idx;
	
	//if (file.is_open())
	//{
		//file << std::scientific;

		//for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		//{
			//m_net.feedforward(dataset.feature_batch(batch));
			//output_batch = m_net.layers.back().states;

			//perturbation.setRandom();
			//perturbation *= 0.001;

			//m_net.feedforward(dataset.feature_batch(batch) + perturbation.matrix());

			//for (int instance = 0; instance < dataset.batch_size; ++instance)
			//{
				//idx = static_cast<long>(instance);

				//for (int i = 0; i < dataset.feature_length(); ++i)
				//{
					//file << std::setw(8);
					//file << m_net.layers[0].states.col(idx)(static_cast<long>(i)) << '\t';
				//}

				//rise = (output_batch.col(idx).matrix() - m_net.layers.back().states.col(idx)).norm();
				//run = perturbation.matrix().col(idx).norm();

				//file << rise / run << '\n';
			//}
		//}
	//}
	//else
	//{
		//open_ascii_escape("red");
		//std::cerr << "\nERROR!: FAILED TO WRITE LYAPUNOV EXPONENT\n";
		//std::cerr << std::endl;
		//close_ascii_escape();
	//}

//}


//template <typename T>
//void Model<T>::save(std::string filename)
//{
	//std::ofstream file(filename.c_str());

	//if (file.is_open())
	//{
		//std::cout << "Saving net to " << filename << '\n';

		//file << "Number of layers:   " << m_net.layers.size() << '\n';

		//for (int i = 0; i < m_net.layers.size(); ++i)
			//file << "Neurons in layer " << i << ": " << m_net.layers[i].biases.rows() << '\n';

		//file << std::scientific;

		//for (int i = 1; i < m_net.layers.size(); ++i)
		//{
			//file << "layer: " << i << '\n';
			//file << "weights: \n";
			//file << m_net.layers[i].weights << '\n';
			//file << "biases: \n";
			//file << m_net.layers[i].biases.col(0) << '\n';
			//file << "--------------\n";
		//}	
	//}
	//else
	//{
		//open_ascii_escape("red");
		//std::cerr << "\nERROR!: FAILED TO SAVE NET\n";
		//std::cerr << std::endl;
		//close_ascii_escape();
	//}
//}

template class Model<double>;
