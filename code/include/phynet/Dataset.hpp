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
#include "Eigen/Dense"
#include "nixio.hpp"
#include "boost/random.hpp"
#include "boost/random/random_device.hpp"

typedef Eigen::Map<const Eigen::MatrixXd> BatchMap;

class Dataset  
{
	public:
		Dataset(std::string input_path, std::size_t ibatch_size, std::string input_type);

		Eigen::MatrixXd feature_batch(std::size_t batch) const;
		Eigen::MatrixXd target_batch(std::size_t batch) const;

		std::size_t feature_length(void) const;
		std::size_t target_length(void) const ;
		std::size_t instances(void) const;
		std::size_t batches(void) const;
		std::size_t batch_size(void) const;
		double energy(std::size_t instance) const;

		Eigen::MatrixXd hamiltonian(std::size_t instance) const;
		Eigen::VectorXd wavefunction(std::size_t instance) const;
		Eigen::VectorXd szdiag(void) const;

		void generate_template_average_file(std::string filename) const;
		void shuffle(void);
		void print_info(void);

		Eigen::VectorXd 
		sparse_ham_times_vec(std::size_t instance, const Eigen::VectorXd &predicted_wavefx) const;
		
	protected:
		void import(std::string input_path);
		void set_topological_structure(void);
		inline int tval(int potential_row_index, int dim) const;

		std::string m_input_type;
		
		std::vector<int> m_header;
		std::vector<int> m_matloc;
		std::vector<double> m_fields;
		std::vector<double> m_szdiag;
		std::vector<double> m_matval;
		std::vector<double> m_wavefx;
		std::vector<double> m_energy;

		std::size_t m_feature_length;
		std::size_t m_target_length;
		std::size_t m_batch_size;

		std::vector<std::size_t> m_indices;
};
	
#endif /* Dataset_hpp */

