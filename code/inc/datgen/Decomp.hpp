///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Decomp.hpp ***                               //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Decomp_hpp
#define Decomp_hpp

#include <vector>
#include <Eigen/Dense>
#include "Instance.hpp"
#include "Operators.hpp"
#include "Fields.hpp"

#define BX_MIN 0.1
#define BX_MAX 2.0
#define J 1.0
#define BZ 0.05

class Decomp  
{
	public:
		Decomp(int num_qubits, int num_transverse, int num_realizations);

		const Eigen::MatrixXd& unitary(int instance) const;
		const Eigen::VectorXd& spectrum(int instance) const;

		void compute(void);

		void prompt_if_file(std::string fpath) const;
		void write(std::string fpath) const;
		void read(std::string fpath);
		void print(int instance) const;

	private:
		int num_qubits, num_transverse, num_realizations;
		std::vector<Instance> instances;
};

	
#endif /* Decomp_hpp */

