///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Gradient.hpp ***                             //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Gradient_hpp
#define Gradient_hpp

#include "Eigen/Dense"

class Gradient  
{
	public:
		Gradient();


	private:
		Eigen::MatrixXd states;
		Eigen::MatrixXd errors;
};

	
#endif /* Gradient_hpp */

