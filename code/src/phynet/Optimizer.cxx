///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Optimizer.hpp ***                             //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Optimizer.hpp>

template <typename T>
Optimizer<T>::Optimizer(std::string optimizer_type, 
		T learning_rate, 
		T decay_rate, 
		T epsilon_conditioner)
	: 
		m_learning_rate(learning_rate), 
		m_decay_rate(decay_rate), 
		m_epsilon_conditioner(epsilon_conditioner)
{
	set_update_pointer(optimizer_type);
}

template <typename T>
void Optimizer<T>::stochastic_gradient_descent(Network<T>& net)
{
	T step_size = m_learning_rate / net.layers[0].states.cols();

	for (std::size_t l = net.layers.size() - 1; l >= 1; --l)
	{
		net.layers[l].weights -= 
			step_size * net.layers[l].errors * net.layers[l-1].states.transpose();
		
		net.layers[l].biases.colwise() -= 
			step_size * net.layers[l].errors.rowwise().sum();
	}	
}

template <typename T>
void Optimizer<T>::adaptive_delta(Network<T>& net)
{
	if (!m_buffers_initialized)
	{
		initialize_accumulation_buffers(net);
		m_buffers_initialized = true;
	}

	compute_gradients(net);
	accumulate_gradient();
	compute_update();
	accumulate_update();
	apply_update(net);
}

template <typename T>
void Optimizer<T>::compute_gradients(const Network<T>& net)
{
	T step_size = m_learning_rate / net.layers[0].states.cols();

	for (std::size_t l = net.layers.size() - 1; l >= 1; --l)
	{
		m_gradients.layers[l].weights = 
			step_size * net.layers[l].errors * net.layers[l-1].states.transpose();

		m_gradients.layers[l].biases.colwise() = 
			step_size * net.layers[l].errors.rowwise().sum();
	}	
}

template <typename T>
void Optimizer<T>::accumulate_gradient(void)
{
	for (std::size_t l = 1; l < m_gradients.layers.size(); ++l)
	{
		m_gradient_accumulator.layers[l].weights = 
			m_decay_rate * m_gradient_accumulator.layers[l].weights.array() +
			(1.0 - m_decay_rate) * m_gradients.layers[l].weights.array().square();
		
		m_gradient_accumulator.layers[l].biases = 
			m_decay_rate * m_gradient_accumulator.layers[l].biases.array() +
			(1.0 - m_decay_rate) * m_gradients.layers[l].biases.array().square();
	}
}

template <typename T>
void Optimizer<T>::compute_update(void)
{
	for (std::size_t l = 1; l < m_gradients.layers.size(); ++l)
	{
		m_updates.layers[l].weights = -
		((m_update_accumulator.layers[l].weights.array() + m_epsilon_conditioner) / 
		(m_gradient_accumulator.layers[l].weights.array() + m_epsilon_conditioner))
		.sqrt() * m_gradients.layers[l].weights.array();

		m_updates.layers[l].biases = -
		((m_update_accumulator.layers[l].biases.array() + m_epsilon_conditioner) / 
		(m_gradient_accumulator.layers[l].biases.array() + m_epsilon_conditioner))
		.sqrt() * m_gradients.layers[l].biases.array();
	}
}

template <typename T>
void Optimizer<T>::accumulate_update(void)
{
	for (std::size_t l = 1; l < m_gradients.layers.size(); ++l)
	{
		m_update_accumulator.layers[l].weights = 
			m_decay_rate * m_update_accumulator.layers[l].weights.array() + 
			(1.0 - m_decay_rate) * m_updates.layers[l].weights.array().square();
		
		m_update_accumulator.layers[l].biases = 
			m_decay_rate * m_update_accumulator.layers[l].biases.array() + 
			(1.0 - m_decay_rate) * m_updates.layers[l].biases.array().square();
	}
}

template <typename T>
void Optimizer<T>::apply_update(Network<T>& net)
{
	for (std::size_t l = 1; l < m_gradients.layers.size(); ++l)
	{
		net.layers[l].weights += m_updates.layers[l].weights;
		net.layers[l].biases += m_updates.layers[l].biases;
	}
}

template <typename T>
void Optimizer<T>::set_update_pointer(std::string optimizer_type)
{
	using namespace std::placeholders;

	if (optimizer_type == "sgd")
	{
		update = std::bind(&Optimizer<T>::stochastic_gradient_descent, this, _1);
	}
	else if (optimizer_type == "adadelta")
	{
		update = std::bind(&Optimizer<T>::adaptive_delta, this, _1);
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "ERROR!: UNSUPPORTED OPTIMIZATION ALGORITHM\n" << std::endl;
		close_ascii_escape();
		exit(-1);
	}

}

template <typename T>
void Optimizer<T>::initialize_accumulation_buffers(const Network<T>& net)
{
	m_updates = net;
	m_gradients = net;
	m_gradient_accumulator = net;
	m_update_accumulator = net;
	
	m_updates.set_zero();
	m_gradients.set_zero();
	m_gradient_accumulator.set_zero();
	m_update_accumulator.set_zero();

	// just to save a little mem and not write new object
	for (std::size_t l = net.layers.size() - 1; l >= 1; --l)
	{
		m_updates.layers[l].errors.resize(0,0);
		m_updates.layers[l].weighted_sum.resize(0,0);

		m_gradients.layers[l].states.resize(0,0);
		m_gradients.layers[l].errors.resize(0,0);
		m_gradients.layers[l].weighted_sum.resize(0,0);
	
		m_gradient_accumulator.layers[l].states.resize(0,0);
		m_gradient_accumulator.layers[l].errors.resize(0,0);
		m_gradient_accumulator.layers[l].weighted_sum.resize(0,0);
		
		m_update_accumulator.layers[l].states.resize(0,0);
		m_update_accumulator.layers[l].errors.resize(0,0);
		m_update_accumulator.layers[l].weighted_sum.resize(0,0);
	}
}



template class Optimizer<float>;
template class Optimizer<double>;
