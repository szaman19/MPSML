///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Layer.cpp ***                               //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Layer.hpp>

template <typename T>
void 
Layer<T>::reserve(int neurons_in_this_layer, int neurons_in_prev_layer, int batch_size)
{
	weights.resize(neurons_in_this_layer, neurons_in_prev_layer);
	biases.resize(neurons_in_this_layer, batch_size);
	states.resize(neurons_in_this_layer, batch_size);
	errors.resize(neurons_in_this_layer, batch_size);
	weighted_sum.resize(neurons_in_this_layer, batch_size);
}

template <typename T>
long Layer<T>::neurons_in_layer(void) const
{
	return states.rows();
}

template <typename T>
void Layer<T>::set_zero(void)
{
	weights.setZero();
	biases.setZero();
	states.setZero();
	errors.setZero();
	weighted_sum.setZero();
}

template <typename T>
void Layer<T>::deep_copy(const Layer<T>& source)
{
	this->weights = source.weights;
	this->biases = source.biases;
	this->states = source.states;
	this->errors = source.errors;
	this->weighted_sum = source.weighted_sum;
}

template <typename T>
Layer<T>& Layer<T>::operator=(const Layer<T>& source)
{
	if (this == &source) return *this;

	deep_copy(source);

	return *this;
}

template <typename T>
void Layer<T>::compute_states(void)
{
	states = weighted_sum.unaryExpr(activation.activation);
}

template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>
Layer<T>::derivative_of_activation_on_weighted_sum(void) const
{
	return weighted_sum.unaryExpr(activation.derivative).array();
}

template class Layer<float>;
template class Layer<double>;

