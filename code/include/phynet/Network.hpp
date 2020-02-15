///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Network.hpp ***                               //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Network_hpp
#define Network_hpp

#include <random>
#include <phynet/Layer.hpp>
#include <phynet/Dataset.hpp>

template <typename T>
class Network
{
	public:
		Network(void) { };
		Network(const Network<T>& source) { deep_copy(source); }
		Network(const Topology<T>& topology);
	
		void init_weights(void);
		void init_biases(T value = 0);

		void validate_topology(const Dataset<T>& dataset) const;
		void set_zero(void);

		Eigen::Matrix<T, 1, Eigen::Dynamic> flat_wandb(void) const;
				
		void feedforward(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& inputs);

		void backpropagate(void);

		Network& operator=(const Network<T>& source);
		
		std::vector<Layer<T>> layers;

	private:
		void deep_copy(const Network<T>& source);
};
	
#endif /* Network_hpp */

