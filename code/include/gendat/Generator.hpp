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

#define BX_MIN 0.02
#define BX_MAX 1.99
#define J 1.0
#define BZ 0.05
#define DISORDER_STRENGTH 0.01

template <typename T>
class Generator  
{
	public:
		Generator(int num_qubits, int num_transverse, int num_realizations);

		void prompt_if_file(std::string fpath) const;
		void set_dump_location(std::string fpath);
		void run(void) const;

	private:
		bool ready_to_dump = false;
		int num_qubits, num_transverse, num_realizations;
		std::string fpath;
};



	
#endif /* Generator_hpp */

