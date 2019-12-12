///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Layer.hpp ***                               //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Layer_hpp
#define Layer_hpp

#include "Eigen/Dense"
#include "Topology.hpp"

class Layer  
{
	public:
		Layer(void) { };
		Layer(Activation activation) : activation(activation) {}

		void compute_states(void);
		void set_zero(void);
		void reserve(std::size_t neurons_in_this_layer, 
		             std::size_t neurons_in_prev_layer,
			         std::size_t batch_size);
			
		long neurons_in_layer(void) const;

		// This needs a better implementation, how not return by value?
		Eigen::ArrayXXd derivative_of_activation_on_weighted_sum(void) const;

		Layer& operator=(const Layer& source);

		Eigen::MatrixXd weights;
		Eigen::MatrixXd biases;
		Eigen::MatrixXd states;
		Eigen::MatrixXd errors;
		Eigen::MatrixXd weighted_sum;

		Activation activation;

	private:
		void deep_copy(const Layer& source);
};

	
#endif /* Layer_hpp */

