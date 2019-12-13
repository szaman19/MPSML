///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** Decomp.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Decomp_hpp
#define Decomp_hpp

#include <vector>
#include <Eigen/Dense>
#include "operators.hpp"

typedef Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes_t;

class Decomp  
{
	public:
		Decomp(int num_qubits, int num_transverse_fields, int num_disorder_realizations);


	//private:
		std::vector<saes_t> wavefx;
		std::vector<Fields> fields;
};

	
#endif /* Decomp_hpp */

