///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Generator.hpp ***                            //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
// modified 2021 by Andrew T Grace, Binghamton CS                            //
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

template <typename T>
class Generator  
{
	public:
		Generator(std::string model);

		void prompt_if_file(std::string fpath) const;

		void run(void) const;
        void runParallel(void) const;

		int dBx, dBz, qubits, replicas;

		double Bx_min, Bx_max, Bz_min, Bz_max, coupling, disorder;

		std::string fpath, model;
};



	
#endif /* Generator_hpp */

