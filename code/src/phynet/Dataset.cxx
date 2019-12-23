///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Dataset.cpp ***                               //
//                                                                           //
// created June 16, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Dataset.hpp>

template <typename T>
Dataset<T>::Dataset(int num_qubits, std::string fpath, int batch_size, int num_instances)
	: batch_size(batch_size), num_qubits(num_qubits), num_instances(num_instances), fpath(fpath)
{ 
	dim = (int)(pow(2, num_qubits) + 0.5);
	init_pos();

	std::cout << "Gathering necessary virtual memory...\n";
	allocate();
	//std::cout << "done\n";

	if (num_instances < batch_size)
	{
		std::cerr << "WARNING: BAD BATCH SIZE...\n";
		exit(-1);
	}

	if ((3*num_qubits + dim + dim*dim)*num_instances*sizeof(T) > 
			boost::filesystem::file_size(fpath))
	{
		std::cerr << "ERROR!: REQUESTED MORE INSTANCES THAN EXIST IN FILE\n";
		exit(-1);
	}

	import(num_instances);
}

template <typename T>
void Dataset<T>::init_pos(void)
{
	std::streamoff offset = (3*num_qubits + dim + dim*dim) * sizeof(T);

	for (int instance = 0; instance < num_instances; ++instance)
		pos.push_back(offset * instance);	
}

template <typename T>
void Dataset<T>::allocate(void)
{
	Batch<T> fields(3 * num_qubits, batch_size);
	Batch<T> energy(dim, batch_size);
	Batch<T> wavefx(dim, batch_size);
	Waves<T> waves;

	int num_training_batches = (int)(num_instances * PERCENT_TRAINING) / batch_size;
	int num_validation_batches = (int)(num_instances * PERCENT_VALIDATION) / batch_size;
	int num_testing_batches = (int)(num_instances * PERCENT_TESTING) / batch_size;

	for (int i = 0; i < dim; ++i) waves.push_back(wavefx);

	for (int i = 0; i < num_training_batches; ++i) 
	{
		training_fields.push_back(fields);
		training_energy.push_back(energy);
		training_wavefx.push_back(waves);
	}

	for (int i = 0; i < num_validation_batches; ++i)
	{
		validation_fields.push_back(fields);
		validation_energy.push_back(energy);
		validation_wavefx.push_back(waves);
	}

	for (int i = 0; i < num_testing_batches; ++i)
	{
		testing_fields.push_back(fields);
		testing_energy.push_back(energy);
		testing_wavefx.push_back(waves);
	}
}

template <typename T>
void Dataset<T>::fill(int pos_lower_idx, int num_batches,
		std::vector<Batch<T>>& fields_var, 
		std::vector<Batch<T>>& energy_var, 
		std::vector<Waves<T>>& wavefx_var)
{
	std::ifstream file(fpath.c_str(), std::ios::binary);  

	if (!file.is_open())
	{
		std::cerr << "ERROR: COULD NOT OPEN " << fpath << '\n';
		exit(-1);
	}

	long fields_bytes = NUM_HAM_PARAM * num_qubits * sizeof(T);
	long wavefx_bytes = dim * sizeof(T);
	long energy_bytes = dim * sizeof(T);

	for (int batch = 0; batch  < num_batches; ++batch)
	{
		for (int instance = 0; instance < batch_size; ++instance)
		{
			// set file pointer 
			file.seekg(pos[pos_lower_idx + instance + batch * batch_size]);
			// pull out fields
			file.read((char*)fields_var[batch].col(instance).data(), fields_bytes);
			// pull out energies
			file.read((char*)energy_var[batch].col(instance).data(), energy_bytes);
			// pull out wavefunctions 	
			for (int vec = 0; vec < dim; ++vec)
				file.read((char*)wavefx_var[batch][vec].col(instance).data(), wavefx_bytes);
		}
	}
}

template <typename T>
void Dataset<T>::import(int num_instances)
{
	std::cout << "Importing " << num_instances << " instances from " << fpath << '\n';
			  
	this->num_instances = num_instances;

	std::random_device rd;	
	std::default_random_engine engine(rd());
	std::shuffle(pos.begin(), pos.end(), engine);

	int lower_training_index = 0;		
	int lower_validation_index = (int)(num_instances * PERCENT_TRAINING);
	int lower_testing_index = lower_validation_index + (int)(num_instances * PERCENT_VALIDATION);

	int num_training_batches = (int)(num_instances * PERCENT_TRAINING) / batch_size;
	int num_validation_batches = (int)(num_instances * PERCENT_VALIDATION) / batch_size;
	int num_testing_batches = (int)(num_instances * PERCENT_TESTING) / batch_size;

	fill(lower_training_index, num_training_batches, 
			training_fields, training_energy, training_wavefx);

	fill(lower_validation_index, num_validation_batches, 
			validation_fields, validation_energy, validation_wavefx);

	fill(lower_testing_index, num_testing_batches, 
			testing_fields, testing_energy, testing_wavefx);
}

template <typename T>
int Dataset<T>::num_training_instances(void) const 
{
	return (int)(num_instances * PERCENT_TRAINING);
}

template <typename T>
int Dataset<T>::num_validation_instances(void) const 
{
	return (int)(num_instances * PERCENT_VALIDATION);
}

template <typename T>
int Dataset<T>::num_testing_instances(void) const 
{
	return (int)(num_instances * PERCENT_TESTING);
}

template <typename T>
int Dataset<T>::feature_length(void) const 
{
	return training_fields[0].rows();
}

template <typename T>
int Dataset<T>::target_length(void) const 
{
	return training_wavefx[0][0].rows();
}

template <typename T>
int Dataset<T>::num_training_batches(void) const 
{
	return training_fields.size();
}

template <typename T>
int Dataset<T>::num_validation_batches(void) const 
{
	return validation_fields.size();
}

template <typename T>
int Dataset<T>::num_eigenvectors(void) const 
{
	return dim;
}

template <typename T>
int Dataset<T>::num_testing_batches(void) const 
{
	return testing_fields.size();
}

template <typename T>
const Batch<T>& Dataset<T>::training_energy_batch(int batch) const 
{
	return training_energy[batch];
}

template <typename T>
const Batch<T>& Dataset<T>::validation_energy_batch(int batch) const 
{
	return training_energy[batch];
}

template <typename T>
const Batch<T>& Dataset<T>::testing_energy_batch(int batch) const 
{
	return training_energy[batch];
}

template <typename T>
const Batch<T>& Dataset<T>::training_feature_batch(int batch) const 
{
	return training_fields[batch];
}
template <typename T>
const Batch<T>& Dataset<T>::validation_feature_batch(int batch) const 
{
	return validation_fields[batch];
}
template <typename T>
const Batch<T>& Dataset<T>::testing_feature_batch(int batch) const 
{
	return testing_fields[batch];
}

template <typename T>
const Batch<T>& Dataset<T>::training_target_batch(int batch, int iwavefx) const 
{
	return training_wavefx[batch][iwavefx];
}

template <typename T>
const Batch<T>& Dataset<T>::validation_target_batch(int batch, int iwavefx) const 
{
	return validation_wavefx[batch][iwavefx];
}

template <typename T>
const Batch<T>& Dataset<T>::testing_target_batch(int batch, int iwavefx) const 
{
	return testing_wavefx[batch][iwavefx];
}

template class Dataset<double>;
