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
#include <phynet/Activation.hpp>

template <typename T>
class Topology  
{
	public:
		Topology(int batch_size);

		int num_layers(void) const;
		int num_neurons_in_layer(int layer) const;
		int batch_size(void) const;

		Activation<T> activation(int layer) const;

		void push_back(int neurons, Activation<T> activations);

	private:
		int m_batch_size;
		std::vector<int> m_neurons_per_layer;	
		std::vector<Activation<T>> m_activations;
};
	
#endif /* Topology_hpp */

