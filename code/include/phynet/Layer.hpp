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

#include <Eigen/Dense>
#include <phynet/Topology.hpp>

template <typename T>
class Layer  
{
	public:
		Layer(void) { };
		Layer(Activation<T> activation) : activation(activation) {}

		void compute_states(void);
		void set_zero(void);

		void reserve(int neurons_in_this_layer, int neurons_in_prev_layer, int batch_size);
			
		long neurons_in_layer(void) const;

		// This needs a better implementation, how not return by value?
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>
		   	derivative_of_activation_on_weighted_sum(void) const;

		Layer<T>& operator=(const Layer<T>& source);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weights;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> biases;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> states;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> errors;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weighted_sum;

		Activation<T> activation;

	private:
		void deep_copy(const Layer<T>& source);
};

	
#endif /* Layer_hpp */

