///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Topology.cpp ***                             //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Topology.hpp"

Topology::Topology(std::size_t batch_size_in)
	:
	m_batch_size(batch_size_in)
{
	// null body
}

std::size_t Topology::num_layers(void) const
{
	return m_neurons_per_layer.size();
}

std::size_t Topology::num_neurons_in_layer(std::size_t layer) const
{
	return m_neurons_per_layer[layer];
}

std::size_t Topology::batch_size(void) const
{
	return m_batch_size;
}

void Topology::push_back(std::size_t neurons, Activation activations)
{
	m_neurons_per_layer.push_back(neurons);
	m_activations.push_back(activations);
}

Activation Topology::activation(std::size_t layer) const
{
	return m_activations[layer];
}
