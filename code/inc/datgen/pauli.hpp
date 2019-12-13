///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** pauli.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef pauli_hpp
#define pauli_hpp

#include <Eigen/Core>

Eigen::Matrix2d eye();

Eigen::Matrix2d sigmax();

Eigen::Matrix2d sigmaz();

Eigen::Matrix2cd sigmay();
	
#endif /* pauli_hpp */

