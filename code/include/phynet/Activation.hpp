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
#include <nixio/nixio.hpp>

#define ELU_ALPHA 0.1

template <typename T>
class Activation  
{
	public:
		Activation(std::string activation_type = "tanh");

		std::function<T(T)> activation;
		std::function<T(T)> derivative;

		static T sigmoid_0(T x);
		static T sigmoid_1(T x);
		static T linear_0(T x);
		static T linear_1(T x);
		static T tanh_0(T x);
		static T tanh_1(T x);
		static T relu_0(T x);
		static T relu_1(T x);
		static T elu_0(T x);
		static T elu_1(T x);
};

	
#endif /* Activation_hpp */

