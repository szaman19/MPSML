///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Network.cpp ***                               //
//                                                                           //
// created June 15, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/Network.hpp>

Network::Network(const Topology &topology)
{
	if (topology.num_layers() < 2)
	{
		open_ascii_escape("yellow");
		std::cerr << "\nWARNING!: NETWORK ONLY HAS ONE LAYER???\n\n";
		close_ascii_escape();
	}
	else
	{
		std::size_t neurons_in_this_layer;
		std::size_t neurons_in_prev_layer;

		for (std::size_t l = 0; l < topology.num_layers(); ++l)
		{
			neurons_in_this_layer = topology.num_neurons_in_layer(l);

			if (l == 0) neurons_in_prev_layer = 1;
			else neurons_in_prev_layer = topology.num_neurons_in_layer(l-1);

			layers.push_back( Layer( topology.activation(l) ) );		

			layers.back().reserve(neurons_in_this_layer, 
					              neurons_in_prev_layer, 
					              topology.batch_size());
		}
	}
	init_weights();
	init_biases();
}

void Network::init_weights()
{
	for (auto& layer : layers) layer.weights.setRandom();
}

void Network::init_biases(double value)
{
	for (auto& layer : layers) layer.biases.fill(value);
}

void Network::validate_topology(const Dataset &dataset) const
{
	if (layers[0].states.rows() != static_cast<long>(dataset.feature_length()))
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: INPUTS INCOMMENSURATE WITH NETWORK TOPOLOGY\n";
		std::cerr << std::endl;
		close_ascii_escape();
		exit(-1);
	}

	if (layers.empty())
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: NETWORK HAS ZERO LAYERS???\n\n";
		close_ascii_escape();
		exit(-1);
	}
	else if (static_cast<long>(dataset.target_length()) 
			 != layers.back().states.rows() ||
			 static_cast<long>(dataset.batch_size())
			 != layers.back().states.cols())
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: TARGETS INCOMMENSURATE WITH NETWORK TOPOLOGY\n\n";
		close_ascii_escape();
		exit(-1);
	}
}

void Network::feedforward(const Eigen::MatrixXd &inputs)
{
	layers[0].states = inputs;

	for (std::size_t layer = 1; layer < layers.size(); ++layer)
	{
		layers[layer].weighted_sum = 
			layers[layer].weights * layers[layer - 1].states + layers[layer].biases;

		layers[layer].compute_states();
	}
}

void Network::backpropagate(void)
{
	for (std::size_t l = layers.size() - 2; l >= 1; --l)
	{
		layers[l].errors = layers[l+1].weights.transpose() * layers[l+1].errors;
		layers[l].errors.array() *= layers[l].derivative_of_activation_on_weighted_sum();
	}	
}

Network& Network::operator=(const Network& source)
{
	if (this == &source) return *this; 

	deep_copy(source);

	return *this;
}

void Network::deep_copy(const Network& source)
{
	this->layers.resize(source.layers.size());

	for (std::size_t i = 0; i < layers.size(); ++i)
		this->layers[i] = source.layers[i];
}

void Network::set_zero(void)
{
	for (std::size_t i = 0; i  < layers.size(); ++i)
	{
		layers[i].set_zero();
	}
}













