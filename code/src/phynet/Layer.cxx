///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Layer.cpp ***                               //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Layer.hpp>

void Layer::reserve(std::size_t neurons_in_this_layer, 
	    		    std::size_t neurons_in_prev_layer,
			   		std::size_t batch_size)
{
	weights.resize(neurons_in_this_layer, neurons_in_prev_layer);
	biases.resize(neurons_in_this_layer, batch_size);
	states.resize(neurons_in_this_layer, batch_size);
	errors.resize(neurons_in_this_layer, batch_size);
	weighted_sum.resize(neurons_in_this_layer, batch_size);
}

long Layer::neurons_in_layer(void) const
{
	return states.rows();
}

void Layer::set_zero(void)
{
	weights.setZero();
	biases.setZero();
	states.setZero();
	errors.setZero();
	weighted_sum.setZero();
}

void Layer::deep_copy(const Layer& source)
{
	this->weights = source.weights;
	this->biases = source.biases;
	this->states = source.states;
	this->errors = source.errors;
	this->weighted_sum = source.weighted_sum;
}

Layer& Layer::operator=(const Layer& source)
{
	if (this == &source) return *this;

	deep_copy(source);

	return *this;
}

void Layer::compute_states(void)
{
	//af::array gpu_weighted_sum(weighted_sum.rows(), weighted_sum.cols(), weighted_sum.data());

	//gpu_weighted_sum = af::tanh(gpu_weighted_sum);

	//typedef Eigen::Map<Eigen::MatrixXd> emap;

	//states = emap(gpu_weighted_sum.host<double>(), weighted_sum.rows(), weighted_sum.cols());

	states = weighted_sum.unaryExpr(activation.activation);
}

Eigen::ArrayXXd Layer::derivative_of_activation_on_weighted_sum(void) const
{
	//af::array gpu_tmp(weighted_sum.rows(), weighted_sum.cols(), weighted_sum.data());

	//gpu_tmp = 1 - af::tanh(gpu_tmp) * af::tanh(gpu_tmp);

	//typedef Eigen::Map<Eigen::ArrayXXd> emap;

	//return emap(gpu_tmp.host<double>(), weighted_sum.rows(), weighted_sum.cols());

	return weighted_sum.unaryExpr(activation.derivative).array();
}

