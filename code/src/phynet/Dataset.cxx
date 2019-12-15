///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Dataset.cpp ***                               //
//                                                                           //
// created June 16, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Dataset.hpp>

Dataset::Dataset(std::string input_path, std::size_t ibatch_size, std::string input_type)
	:
	m_input_type(input_type), 
	m_batch_size(ibatch_size)
{
	import(input_path);	
	m_target_length = static_cast<std::size_t>(m_header[2]);
	set_topological_structure();
	m_indices.resize(m_header[0]);
	std::iota(m_indices.begin(), m_indices.end(), 0);
}

void Dataset::import(std::string input_path)
{
	std::string header_path = input_path + "header.bin";
	std::string szdiag_path = input_path + "szdiag.bin";
	std::string fields_path = input_path + "fields.bin";
	std::string matloc_path = input_path + "matloc.bin";
	std::string matval_path = input_path + "matval.bin";
	std::string wavefx_path = input_path + "wavefx.bin";
	std::string energy_path = input_path + "energy.bin";

	std::ifstream header_file(header_path.c_str(), std::ios::binary);
	std::ifstream szdiag_file(szdiag_path.c_str(), std::ios::binary);
	std::ifstream fields_file(fields_path.c_str(), std::ios::binary);
	std::ifstream matloc_file(matloc_path.c_str(), std::ios::binary);
	std::ifstream matval_file(matval_path.c_str(), std::ios::binary);
	std::ifstream wavefx_file(wavefx_path.c_str(), std::ios::binary);
	std::ifstream energy_file(energy_path.c_str(), std::ios::binary);

	std::cout << "Reading dataset from " << input_path << '\n';

	if (!header_file.is_open() || 
		!szdiag_file.is_open() ||
		!fields_file.is_open() ||
		!matloc_file.is_open() || 
		!matval_file.is_open() || 
		!wavefx_file.is_open() || 
		!energy_file.is_open())
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO IMPORT DATA FROM " 
			<< input_path << '\n' << std::endl;
		close_ascii_escape();
		exit(-1);
	}

	m_header.reserve(3);
	header_file.read(reinterpret_cast<char*>(m_header.data()), static_cast<long>(3*sizeof(int)));

	std::size_t tmp = static_cast<std::size_t>(m_header[0]);
	std::size_t nnz = static_cast<std::size_t>(m_header[1]);
	std::size_t dim = static_cast<std::size_t>(m_header[2]);

	m_fields.resize(tmp);
	m_szdiag.resize(dim);
	m_matloc.resize(tmp * nnz);
	m_matval.resize(tmp * nnz);
	m_wavefx.resize(tmp * dim);
	m_energy.resize(tmp);

	long fields_bytes = static_cast<long>(tmp * sizeof(double));
	long szdiag_bytes = static_cast<long>(dim * sizeof(double));
	long matloc_bytes = static_cast<long>(tmp * nnz * sizeof(int));
	long matval_bytes = static_cast<long>(tmp * nnz * sizeof(double));
	long wavefx_bytes = static_cast<long>(tmp * dim * sizeof(double));
	long energy_bytes = static_cast<long>(tmp * sizeof(double));
	
	fields_file.read(reinterpret_cast<char*>(m_fields.data()), fields_bytes);
	szdiag_file.read(reinterpret_cast<char*>(m_szdiag.data()), szdiag_bytes);
	matloc_file.read(reinterpret_cast<char*>(m_matloc.data()), matloc_bytes);
	matval_file.read(reinterpret_cast<char*>(m_matval.data()), matval_bytes);
	wavefx_file.read(reinterpret_cast<char*>(m_wavefx.data()), wavefx_bytes);
	energy_file.read(reinterpret_cast<char*>(m_energy.data()), energy_bytes);
}

std::size_t Dataset::feature_length(void) const
{
	return m_feature_length;
}

std::size_t Dataset::target_length(void) const
{
	return m_target_length;
}

std::size_t Dataset::instances(void) const
{
	return static_cast<std::size_t>(m_header[0]);
}

std::size_t Dataset::batches(void) const
{
	return static_cast<std::size_t>(m_header[0])/m_batch_size;
}

std::size_t Dataset::batch_size(void) const
{
	return m_batch_size;
}

Eigen::VectorXd Dataset::szdiag(void) const
{
	return Eigen::Map<const Eigen::VectorXd>(m_szdiag.data(), m_header[2]);
}

Eigen::MatrixXd Dataset::feature_batch(std::size_t batch) const
{
	std::size_t idx = batch * m_feature_length * m_batch_size;

	return BatchMap(&m_fields[idx], 
			static_cast<long>(m_feature_length), 
			static_cast<long>(m_batch_size));
}

Eigen::MatrixXd Dataset::target_batch(std::size_t batch) const
{
	std::size_t idx = batch * m_target_length * m_batch_size;

	return BatchMap(&m_wavefx[idx],
			static_cast<long>(m_target_length),
			static_cast<long>(m_batch_size));
}

void Dataset::shuffle(void)
{
	// will ensure different ordering everytime
	boost::random::random_device rd;
	boost::random::mt19937 engine(rd);

	// shuffle the vector of instance labels 
	std::shuffle(m_indices.begin(), m_indices.end(), engine);

	// right now have very bad just brute force shuffle algorithm
	std::vector<double> tmp(m_fields);

	for (std::size_t i = 0; i < instances(); ++i) 
		tmp[i] = m_fields[ m_indices[i] ];
	
	m_fields = tmp;

	for (std::size_t i = 0; i < instances(); ++i)
		tmp[i] = m_energy[ m_indices[i] ];

	m_energy = tmp;

	std::size_t nnz = static_cast<std::size_t>(m_header[1]);
	std::size_t dim = static_cast<std::size_t>(m_header[2]);

	tmp.resize(instances() * nnz);

	for (std::size_t i = 0; i < instances(); ++i)
		for (std::size_t j = 0; j < nnz; ++j)
			tmp[i * nnz + j] = m_matval[ m_indices[i] * nnz + j];

	m_matval = tmp;
	
	tmp.resize(instances() * dim);

	for (std::size_t i = 0; i < instances(); ++i)
		for (std::size_t j = 0; j < dim; ++j)
			tmp[i * dim + j] = m_wavefx[ m_indices[i] * dim + j];

	std::vector<int> tmp2(m_matloc);
	
	for (std::size_t i = 0; i < instances(); ++i)
		for (std::size_t j = 0; j < nnz; ++j)
			tmp2[i * nnz + j] = m_matloc[ m_indices[i] * nnz + j];

	m_matloc = tmp2;
}

void Dataset::set_topological_structure(void)
{
	if (m_input_type == "fields") 
		m_feature_length = m_fields.size()/instances();
	else if (m_input_type == "hamiltonian_nz")
		m_feature_length = static_cast<std::size_t>(m_header[1]);
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: INPUT TYPE " 
			<< m_input_type << " INVALID\n" << std::endl;
		close_ascii_escape();
		exit(-1);
	}
}

Eigen::MatrixXd Dataset::hamiltonian(std::size_t instance) const
{
	std::size_t nnz = static_cast<std::size_t>(m_header[1]);
	std::size_t dim = static_cast<std::size_t>(m_header[2]);
	std::size_t lower_bound = instance * nnz;	

	int d = m_header[2];

	Eigen::MatrixXd ham = Eigen::MatrixXd::Zero(d, d);
	std::vector<double> uta;
	uta.resize(dim * (dim + 1) / 2);

	std::size_t index;
	for (std::size_t i = 0; i < nnz; ++i)
	{
		index = lower_bound + i;
		uta[static_cast<std::size_t>(m_matloc[index])] = m_matval[index];
	}

	index = 0;
	for (int i = 0; i < d; ++i)
	{
		for (int j = i; j < d; ++j)
		{
			ham(i,j) = uta[index];
			ham(j,i) = uta[index++];
		}
	}
	return ham;
}

void Dataset::generate_template_average_file(std::string filename) const 
{
	std::ofstream file(filename);

	int dim = m_header[2];

	double avg_of_value, avg_of_square;
	double coeff_squared, coeff;

	Eigen::VectorXd szdiag = this->szdiag();
	Eigen::VectorXd wavefx;

	if (file.is_open())
	{
		file << "#Field   <Sz>   <Sz^2>   <Sz>^2 \n";
		file << std::scientific;

		for (std::size_t instance = 0; instance < instances(); ++instance)
		{
			file << std::setw(8);
			file << m_fields[instance] << '\t';

			wavefx = this->wavefunction(instance);
			avg_of_value = 0;
			avg_of_square = 0;

			for (int i = 0; i < dim; ++i)
			{
				coeff = wavefx(i);
				coeff_squared = coeff * coeff;
				avg_of_value += szdiag(i) * coeff_squared; 
				avg_of_square += szdiag(i) * szdiag(i) * coeff_squared;
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
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: FAILED TO OPEN " << filename << '\n';
		std::cerr << std::endl;
		close_ascii_escape();
		exit(-1);
	}
}

Eigen::VectorXd Dataset::wavefunction(std::size_t instance) const
{
	std::size_t dim = static_cast<std::size_t>(m_header[2]);
	std::size_t lower_bound = instance * dim;	

	int d = m_header[2];

	return Eigen::Map<const Eigen::VectorXd>(&m_wavefx[lower_bound], d); 
}

inline int Dataset::tval(int pot_i, int dim) const
{
	return (pot_i + 1) * dim - ((pot_i*pot_i + pot_i + 2)/2);
}

void Dataset::print_info(void)
{
	std::cout << "INT: " << m_header[0] << '\n';
	std::cout << "NNZ: " << m_header[1] << '\n';
	std::cout << "DIM: " << m_header[2] << '\n';
	std::cout << "BAT: " << m_batch_size << '\n';
	std::cout << "TYP: " << m_input_type << '\n';
	std::cout << "FET: " << m_feature_length << '\n';
	std::cout << "TAR: " << m_target_length << '\n';
	std::cout << "FEL: " << m_fields.size() << '\n';
}

Eigen::VectorXd 
Dataset::sparse_ham_times_vec(std::size_t instance, const Eigen::VectorXd &predicted_wavefx) const
{
	int nnz = m_header[1];
	int dim = m_header[2];

	Eigen::VectorXd resultant = Eigen::VectorXd::Zero(dim);

	int pre_row_ind = 0;
	int act_row_ind, act_col_ind;

	std::size_t offset = instance * static_cast<std::size_t>(nnz);

	for (std::size_t nz = 0; nz < static_cast<std::size_t>(nnz); ++nz)
	{
		for (int pot_row_ind = pre_row_ind; pot_row_ind < dim; ++pot_row_ind)
		{
			if (m_matloc[nz + offset] <= tval(pre_row_ind, dim))
			{
				act_row_ind = pot_row_ind;
				pre_row_ind = pot_row_ind;
				act_col_ind = m_matloc[nz + offset] - (tval(pot_row_ind - 1, dim) + 1) + act_row_ind;

				resultant[act_row_ind] += m_matval[nz + offset] * predicted_wavefx[act_col_ind];
				if (act_row_ind != act_col_ind)
					resultant[act_col_ind] += m_matval[nz + offset] * predicted_wavefx[act_row_ind];

				break;
			}
		}
	}

	return resultant;
}

double Dataset::energy(std::size_t instance) const 
{
	return m_energy[instance];
}
