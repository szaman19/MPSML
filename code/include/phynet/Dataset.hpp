///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Dataset.hpp ***                               //
//                                                                           //
// created June 16, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Dataset_hpp
#define Dataset_hpp

#include <vector>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include <nixio/nixio.hpp>
#include <gendat/Fields.hpp>
#include <boost/filesystem.hpp>
#include <gendat/Operators.hpp>

#define PERCENT_TRAINING 0.01
#define PERCENT_VALIDATION 0.1
#define PERCENT_TESTING 0.90
#define NUM_HAM_PARAM 3

template <typename T>
using Batch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using Waves = std::vector<Batch<T>>;

template <typename T>
class Dataset  
{
	public:
		Dataset(int num_qubits, std::string fpath, int batch_size, 
				std::string input, const Operators<T>& operators);

		const Batch<T>& training_energy_batch(int batch) const;
		const Batch<T>& validation_energy_batch(int batch) const;
		const Batch<T>& testing_energy_batch(int batch) const;

		const Batch<T>& training_feature_batch(int batch) const;
		const Batch<T>& validation_feature_batch(int batch) const;
		const Batch<T>& testing_feature_batch(int batch) const;

		const Batch<T>& training_target_batch(int batch, int iwavefx) const;
		const Batch<T>& validation_target_batch(int batch, int iwavefx) const;
		const Batch<T>& testing_target_batch(int batch, int iwavefx) const;

		int feature_length(void) const;
		int target_length(void) const ;

		int num_training_batches(void) const;
		int num_validation_batches(void) const;
		int num_testing_batches(void) const;
		int num_eigenvectors(void) const;

		int num_training_instances(void) const;
		int num_validation_instances(void) const;
		int num_testing_instances(void) const;

		void import(void);

		void shuffle(void);

		const int batch_size;
		
	private:
		void allocate(void);
		void init_pos(void);

		void fill_hamnze(const std::vector<Batch<T>>& fields_var, 
				std::vector<Batch<T>>& hamnze_var);

		void fill(int pos_lower_idx, int num_batches,
				std::vector<Batch<T>>& fields_var, 
				std::vector<Batch<T>>& energy_var, 
				std::vector<Waves<T>>& wavefx_var);

		int num_qubits, num_instances, dim;

		std::string fpath, input;

		Operators<T> operators;

		std::vector<std::streampos> pos;

		std::vector<Batch<T>> training_fields;
		std::vector<Batch<T>> training_hamnze;
		std::vector<Batch<T>> training_energy;
		std::vector<Waves<T>> training_wavefx;

		std::vector<Batch<T>> validation_fields;
		std::vector<Batch<T>> validation_hamnze;
		std::vector<Batch<T>> validation_energy;
		std::vector<Waves<T>> validation_wavefx;

		std::vector<Batch<T>> testing_fields;
		std::vector<Batch<T>> testing_hamnze;
		std::vector<Batch<T>> testing_energy;
		std::vector<Waves<T>> testing_wavefx;
};





#endif /* Dataset_hpp */

