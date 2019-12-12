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
#include "Layer.hpp"
#include "Dataset.hpp"

class Network
{
	public:
		Network(void) { };
		Network(const Network& source) { deep_copy(source); }
		Network(const Topology &topology);
	
		void init_weights(void);
		void init_biases(double value = 0);

		void validate_topology(const Dataset &dataset) const;
		void set_zero(void);
				
		void feedforward(const Eigen::MatrixXd &inputs);
		void backpropagate(void);

		Network& operator=(const Network& source);
		
		std::vector<Layer> layers;

	private:
		void deep_copy(const Network& source);
};

	
#endif /* Network_hpp */

