///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Generator.hpp ***                            //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Generator_hpp
#define Generator_hpp

#include <vector>
#include <Eigen/Dense>
#include "Instance.hpp"
#include "Operators.hpp"
#include "Fields.hpp"
#include <chrono>

#define BX_MIN 0.00
#define BX_MAX 2.00
#define J 1.0
#define BZ 0.01
#define DISORDER_STRENGTH 0.01

template <typename T>
class Generator  
{
	public:
		Generator(std::string model, int num_qubits, int num_transverse, int num_realizations);

		void prompt_if_file(std::string fpath) const;
		void set_dump_location(std::string fpath);
		void run(void) const;

	private:
		bool ready_to_dump = false;
		int num_qubits, num_transverse, num_realizations;
		std::string fpath, model;
};



	
#endif /* Generator_hpp */

