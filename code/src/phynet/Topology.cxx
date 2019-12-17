///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Topology.cpp ***                             //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Topology.hpp>

template <typename T>
Topology<T>::Topology(int batch_size_in)
	:
	m_batch_size(batch_size_in)
{
	// null body
}

template <typename T>
int Topology<T>::num_layers(void) const
{
	return m_neurons_per_layer.size();
}

template <typename T>
int Topology<T>::num_neurons_in_layer(int layer) const
{
	return m_neurons_per_layer[layer];
}

template <typename T>
int Topology<T>::batch_size(void) const
{
	return m_batch_size;
}

template <typename T>
void Topology<T>::push_back(int neurons, Activation<T> activations)
{
	m_neurons_per_layer.push_back(neurons);
	m_activations.push_back(activations);
}

template <typename T>
Activation<T> Topology<T>::activation(int layer) const
{
	return m_activations[layer];
}

template class Topology<float>;
template class Topology<double>;
