///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Topology.hpp ***                              //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Topology_hpp
#define Topology_hpp

#include <vector>
#include "Activation.hpp"

class Topology  
{
	public:
		Topology(std::size_t batch_size);

		std::size_t num_layers(void) const;
		std::size_t num_neurons_in_layer(std::size_t layer) const;
		std::size_t batch_size(void) const;

		Activation activation(std::size_t layer) const;

		void push_back(std::size_t neurons, Activation activations);

	private:
		std::size_t m_batch_size;
		std::vector<std::size_t> m_neurons_per_layer;	
		std::vector<Activation> m_activations;
};
	
#endif /* Topology_hpp */

