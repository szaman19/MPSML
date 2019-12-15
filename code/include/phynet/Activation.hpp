///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Activation.hpp ***                             //
//                                                                           //
// created June 17, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Activation_hpp
#define Activation_hpp

#include <cmath>
#include <string>
#include <functional>
#include "nixio.hpp"

#define ELU_ALPHA 0.1

class Activation  
{
	public:
		Activation(std::string activation_type = "tanh");

		std::function<double(double)> activation;
		std::function<double(double)> derivative;

		static double sigmoid_0(double x);
		static double sigmoid_1(double x);
		static double linear_0(double x);
		static double linear_1(double x);
		static double tanh_0(double x);
		static double tanh_1(double x);
		static double relu_0(double x);
		static double relu_1(double x);
		static double elu_0(double x);
		static double elu_1(double x);
};

	
#endif /* Activation_hpp */

