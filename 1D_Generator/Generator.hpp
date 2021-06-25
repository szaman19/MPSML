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

template <typename T>
class Generator  
{
	public:
		Generator();

		void run(void) const;

		int dBx, dBz, qubits;

		double Bx_min, Bx_max, Bz_min, Bz_max, coupling;

		std::string fpath;
};



	
#endif /* Generator_hpp */

