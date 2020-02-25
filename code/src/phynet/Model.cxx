///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Model.hpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Model.hpp>
#include <boost/circular_buffer.hpp>

template <typename T>
Model<T>::
Model(const std::vector<Network<T>>& networks, const Loss<T>& loss, const Optimizer<T>& optimizer) 
	: networks(networks), loss(loss), optimizer(optimizer) 
{ 
 	// Null constructor body 
}

template <typename T>
void Model<T>::append_wandb_for_radviz(std::string fpath)
{
	std::ofstream file(fpath.c_str(), std::ios::app | std::ios::binary);

	Eigen::Matrix<T, 1, Eigen::Dynamic> out = this->networks[0].flat_wandb();

	file.write((char*)out.data(), out.size() * sizeof(T));
}

template <typename T>
void Model<T>::write_entanglement_entropy(const Dataset<T>& dataset, std::string fpath)
{
	std::ofstream file(fpath);

	if (file.is_open())
	{
		Eigen::Matrix<T, Eigen::Dynamic, 1> psi, tar;
		int n = dataset.operators.num_qubits;
		file << std::scientific;

		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		{
			for (std::size_t net = 0; net < this->networks.size(); ++net)
				this->networks[net].feedforward(dataset.testing_feature_batch(batch));

			for (int instance = 0; instance < dataset.batch_size; ++instance)
			{
				file << dataset.testing_feature_batch(batch).col(instance)(n) << '\t';
				file << dataset.testing_feature_batch(batch).col(instance)(2*n) << '\t';
				for (std::size_t net = 0; net < this->networks.size(); ++net)
				{
					psi = this->networks[net].layers.back().states.col(instance);
					psi /= psi.norm();
					tar = dataset.testing_target_batch(batch, net).col(instance);
					file << entanglement_entropy(tar) << '\t';
					file << entanglement_entropy(psi) << '\t';
				}
				file << '\n';
			}
		}
	}
	else
	{
		std::cerr << "Could not open " << fpath << " for writing\n";
		exit(-1);
	}
}

template <typename T>
T Model<T>::entanglement_entropy(Eigen::Matrix<T, Eigen::Dynamic, 1> psi) const
{
	T out = 0;
	
	int full_dim = psi.size();
	int num_qubits = (int)std::log2(full_dim);
	int redu_dim = std::pow(2, num_qubits/2);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		redu_den(redu_dim, redu_dim), full_den(full_dim, full_dim);

	full_den = psi * psi.transpose();
	redu_den.setZero();

	for (int k = 0; k < redu_dim; ++k)
	{
		for (int l = 0; l < redu_dim; ++l)
		{
			for (int j = 0; j < redu_dim; ++j)
			{
				int row = k*redu_dim + j;
				int col = l*redu_dim + j;
				redu_den(k,l) += full_den(row, col);
			}
		}
	}

	Eigen::SelfAdjointEigenSolver< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > saes(redu_den);

	if (saes.info() != Eigen::Success)
	{
		open_ascii_escape("red");
		std::cout << "ERROR: DIAGONALIZATION OF DENSITY MATRIX FAILED" << std::endl;
		close_ascii_escape();
		exit(-1);
	}

	for (int i = 0; i < redu_dim; ++i)
	{
		if (saes.eigenvalues()(i) - 1e-16 < 0) 
		{
			//open_ascii_escape("yellow");
			//std::cout << "WARNING: POSSIBLE DOMAIN ERROR IN LOGARITHM" << std::endl;
			//close_ascii_escape();
			continue;
		}
		else
		{
			out += saes.eigenvalues()(i) * std::log2(saes.eigenvalues()(i));
		}
	}

	return -out;
}

template <typename T>
void Model<T>::write_radial_visualization(const Dataset<T> &dataset, std::string fpath)
{
	std::ofstream file(fpath);

	if (file.is_open())
	{
		file << std::scientific;

		int row = 0;
		int dim = dataset.num_eigenvectors();
		double c = 2 * 3.14159 / dim;
		double radius = 1;

		auto anchor = [c](int j, double r) -> Eigen::RowVector2d
		{
			return Eigen::RowVector2d(r*cos((j-0)*c), r*sin((j-0)*c));
		};

		std::ofstream ank_file("anchors.dat");
		ank_file << std::scientific;
		for (int j = 0; j < dim; ++j) ank_file << anchor(j, radius) << '\n';

		Eigen::RowVector2d out_true = Eigen::Vector2d::Zero();
		Eigen::RowVector2d out_pred = Eigen::Vector2d::Zero();
		Eigen::Matrix<T, Eigen::Dynamic, 1> psi, tar;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> x(dataset.num_testing_instances(), dim);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> w(dataset.num_testing_instances(), dim);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y(dataset.num_testing_instances(), dim);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> z(dataset.num_testing_instances(), dim);

		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		{
			for (std::size_t net = 0; net < this->networks.size(); ++net)
				this->networks[net].feedforward(dataset.testing_feature_batch(batch));

			for (int instance = 0; instance < dataset.batch_size; ++instance)
			{
				for (std::size_t net = 0; net < this->networks.size(); ++net)
				{
					psi = this->networks[net].layers.back().states.col(instance);
					y.row(row) = psi/psi.norm();
					x.row(row++) = dataset.testing_target_batch(batch, 0).col(instance);
				}
			}
		}

		double max, min;

		max = w.maxCoeff();
		min = w.minCoeff();

		for (int i = 0; i < w.cols(); ++i)
		{
			for (int j = 0; j < w.rows(); ++j)
			{
				w(j,i) = ( w(j,i) - min ) / (max - min);
			}
		}

		max = y.maxCoeff();
		min = y.minCoeff();

		for (int i = 0; i < w.cols(); ++i)
		{
			for (int j = 0; j < w.rows(); ++j)
			{
				y(j,i) = ( y(j,i) - min ) / (max - min);
			}
		}

		for (int i = 0; i < w.rows(); ++i)
		{
			for (int j = 0; j < w.cols(); ++j)
			{
				w(i,j) = x(i,j) / x.row(i).sum();
				z(i,j) = y(i,j) / y.row(i).sum();
			}
		}

		for (int i = 0; i < x.rows(); ++i)
		{
			for (int j = 0; j < x.cols(); ++j)
			{
				out_true += w(i, j) * anchor(j, radius);
				out_pred += z(i, j) * anchor(j, radius);
			}	

			file << out_true << '\t' << out_pred << '\n';
			out_true.setZero();
			out_pred.setZero();
		}
		
		//std::cout << "wrote radviz to " << fpath << '\n';

	}
	else
	{
		std::cerr << "Failed to open " << fpath << " for write" << std::endl;
		exit(-1);
	}

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

	for (std::size_t net = 0; net < networks.size(); ++net)
	{
		for (int batch = 0; batch < dataset.num_validation_batches(); ++batch)
		{
			networks[net].feedforward(dataset.validation_feature_batch(batch));

			tmp += (networks[net].layers.back().states 
					- dataset.validation_target_batch(batch, net)).array().square();
		}
	}
	
	return tmp.mean() / (networks.size() * dataset.num_validation_batches());	
}

template <typename T>
void Model<T>::print_inference_time(const Dataset<T> &dataset)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

	start = std::chrono::high_resolution_clock::now();

	for (std::size_t net = 0; net < networks.size(); ++net)
		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
			networks[net].feedforward(dataset.testing_feature_batch(batch));

	stop = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = stop - start;
	std::cout << "Inference time: " << elapsed.count() << " seconds" << std::endl;
}

template <typename T>
void Model<T>::print_average_sz_error(const Dataset<T> &dataset)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sz = dataset.operators.magnetization();
	Eigen::Matrix<T, Eigen::Dynamic, 1> psi, tar;

	double out = 0;
	double t1, t2;

	for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
	{
		for (std::size_t net = 0; net < this->networks.size(); ++net)
			this->networks[net].feedforward(dataset.testing_feature_batch(batch));

		for (int instance = 0; instance < dataset.batch_size; ++instance)
		{
			for (std::size_t net = 0; net < this->networks.size(); ++net)
			{
				psi = this->networks[net].layers.back().states.col(instance);
				psi /= psi.norm();
				tar = dataset.testing_target_batch(batch, net).col(instance);
				t1 = psi.transpose() * Sz * psi;
				t2 = tar.transpose() * Sz * tar;
				out += std::fabs(t2 - t1);
			}
		}
	}

	std::cout << "Average Sz error: " 
		<< out / (dataset.batch_size * dataset.num_testing_batches() * networks.size())
		<< std::endl;
}

template <typename T>
void Model<T>::write_magnetization(const Dataset<T> &dataset, std::string fpath)
{
	std::ofstream file(fpath);

	if (file.is_open())
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sz = dataset.operators.magnetization();
		Eigen::Matrix<T, Eigen::Dynamic, 1> psi, tar;
		int n = dataset.operators.num_qubits;
		file << std::scientific;

		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		{
			for (std::size_t net = 0; net < this->networks.size(); ++net)
				this->networks[net].feedforward(dataset.testing_feature_batch(batch));

			for (int instance = 0; instance < dataset.batch_size; ++instance)
			{
				file << dataset.testing_feature_batch(batch).col(instance)(n) << '\t';
				file << dataset.testing_feature_batch(batch).col(instance)(2*n) << '\t';
				for (std::size_t net = 0; net < this->networks.size(); ++net)
				{
					psi = this->networks[net].layers.back().states.col(instance);
					tar = dataset.testing_target_batch(batch, net).col(instance);
					psi /= psi.norm();
					file << tar.transpose() * Sz * tar << '\t';
					file << psi.transpose() * Sz * psi << '\t';
				}
				file << '\n';	
			}
		}

		//std::cout << "wrote magnetization to " << fpath << '\n';
	}
	else
	{
		std::cerr << "Could not open " << fpath << " for writing\n";
		exit(-1);
	}
}

template <typename T>
void Model<T>::print_average_overlap(const Dataset<T>& dataset)
{
	int dim = dataset.num_eigenvectors();

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_preds(dim, dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_targs(dim, dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> normalized_final(dim, dim);

	normalized_final.setZero();

	for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
	{
		for (std::size_t net = 0; net < networks.size(); ++net)
			networks[net].feedforward(dataset.testing_feature_batch(batch));

		for (int instance = 0; instance < dataset.batch_size; ++instance)
		{
			for (std::size_t net = 0; net < networks.size(); ++net)
			{
				normalized_targs.col(net) = 
					dataset.testing_target_batch(batch, net).col(instance);

				normalized_preds.col(net) = 
					networks[net].layers.back().states.col(instance) / 
					networks[net].layers.back().states.col(instance).norm();
			}

			normalized_final += normalized_preds.transpose() * normalized_targs;
		}
	}

	normalized_final /= dataset.num_testing_instances();

	if (networks.size() == 1)
	{
		std::cout << "Average overlap: ";
		std::cout << normalized_final(0,0) << '\n';
	}
	else if ((int)networks.size() == dim && dim <= 16)
	{
		std::cout << "\nAverage overlap: \n";
		pretty_print(normalized_final);
	}
	else
	{
		std::cout << "Mean: " << normalized_final.mean() << '\n';
		std::cout << "Determinant: " << normalized_final.determinant() << '\n';
		std::cout << "Norm: " << (1.0 / std::sqrt(dim)) * normalized_final.norm() << "\n\n";
	}
}

template <typename T>
void Model<T>::write_overlap(const Dataset<T>& dataset, std::string fpath)
{
	std::ofstream file(fpath);

	if (file.is_open())
	{
		int n = dataset.operators.num_qubits;
		Eigen::Matrix<T, Eigen::Dynamic, 1> psi, tar;
		file << std::scientific;

		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		{
			for (std::size_t net = 0; net < networks.size(); ++net)
				networks[net].feedforward(dataset.testing_feature_batch(batch));

			for (int instance = 0; instance < dataset.batch_size; ++instance)
			{
				file << dataset.testing_feature_batch(batch).col(instance)(n) << '\t';
				file << dataset.testing_feature_batch(batch).col(instance)(2*n) << '\t';
				for (std::size_t net = 0; net < this->networks.size(); ++net)
				{
					psi = this->networks[net].layers.back().states.col(instance);
					tar = dataset.testing_target_batch(batch, net).col(instance);
					psi /= psi.norm();
					file << psi.transpose() * tar << '\t';
				}
				file << '\n';	
			}
		}
		//std::cout << "wrote overlaps to " << fpath << '\n';
	}
	else
	{
		std::cerr << "Could not open " << fpath << " for writing\n";
		exit(-10);
	}
}

template <typename T>
void Model<T>::write_lyapunov_estimate(const Dataset<T> &dataset, std::string filename, int epoch)
{
	std::ofstream file(filename + "-epoch=" + std::to_string(epoch) + ".dat");

	if (file.is_open())
	{
		file << std::scientific;

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			perturbation(dataset.feature_length(), dataset.batch_size);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			output_batch(dataset.target_length(), dataset.batch_size);

		T rise, run;
		int n = dataset.operators.num_qubits;

		for (int batch = 0; batch < dataset.num_testing_batches(); ++batch)
		{
			for (std::size_t net = 0; net < networks.size(); ++net)
				networks[net].feedforward(dataset.testing_feature_batch(batch));

			output_batch = networks[0].layers.back().states;

			perturbation.setRandom();
			perturbation *= 0.0001;
			
			for (std::size_t net = 0; net < networks.size(); ++net)
				networks[net].feedforward(dataset.testing_feature_batch(batch) + perturbation);

			for (int instance = 0; instance < dataset.batch_size; ++instance)
			{
				rise = (output_batch.col(instance) - 
						networks[0].layers.back().states.col(instance)).norm();

				run = perturbation.col(instance).norm();

				file << dataset.testing_feature_batch(batch).col(instance)(n) << '\t';
				file << dataset.testing_feature_batch(batch).col(instance)(2*n) << '\t';
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
