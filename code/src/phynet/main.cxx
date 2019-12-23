///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** main.cpp ***                               //
//                                                                           //
// created June 13, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <limits>
#include <string>
#include <boost/circular_buffer.hpp>
#include <phynet/Model.hpp>
#include <phynet/Parser.hpp>
#include <phynet/Dataset.hpp>
#include <gendat/Operators.hpp>

int main( int argc, char *argv[] )
{
	std::string param_file;

	if (argc == 2) param_file = argv[1];
	else param_file = "etc/phynet.in";

	Parser parser(param_file);
	
	const std::string loss_type = 
		parser.value_of_key("loss_type").empty() ? 
		"bb" : parser.value_of_key("loss_type");
	
	const int batch_size = 
		parser.value_of_key("batch_size").empty() ? 
		100 : std::atoi(parser.value_of_key("batch_size").c_str());

	const std::string cutout_option = 
		parser.value_of_key("cutout_option").empty() ? 
		"true" : parser.value_of_key("cutout_option");
	
	const int epochs = 
		parser.value_of_key("epochs").empty() ? 
		1000 : std::atoi(parser.value_of_key("epochs").c_str());
	
	const std::string data_root = 
		parser.value_of_key("data_root").empty() ? 
		"data" : parser.value_of_key("data_root");

	const std::string hidden_activation_type = 
		parser.value_of_key("hidden_activation").empty() ? 
		"tanh" : parser.value_of_key("hidden_activation");
	
	const int hidden_layer_size =
		parser.value_of_key("hidden_layer_size").empty() ? 
		100 : std::atoi(parser.value_of_key("hidden_layer_size").c_str());
	
	const double inactive_threshold = 
		parser.value_of_key("inactive_threshold").empty() ?
		1e-6 : std::atof(parser.value_of_key("inactive_threshold").c_str());

	const double lagrange_multiplier = 
		parser.value_of_key("lagrange_multiplier").empty() ?
		0.1 : std::atof(parser.value_of_key("lagrange_multiplier").c_str());
	
	const double random_domain_bound = 
		parser.value_of_key("random_domain_bound").empty() ?
		0.1 : std::atof(parser.value_of_key("random_domain_bound").c_str());

	const double trade_off_parameter = 
		parser.value_of_key("trade_off_parameter").empty() ?
		0.1 : std::atof(parser.value_of_key("trade_off_parameter").c_str());
	
	const double decay_rate = 
		parser.value_of_key("decay_rate").empty() ?
		0.9 : std::atof(parser.value_of_key("decay_rate").c_str());
	
	const double epsilon_conditioner = 
		parser.value_of_key("epsilon_conditioner").empty() ?
		0.95 : std::atof(parser.value_of_key("epsilon_conditioner").c_str());

	const double learning_rate = 
		parser.value_of_key("learning_rate").empty() ?
		0.1 : std::atof(parser.value_of_key("learning_rate").c_str());
	
	const std::string optimizer_type = 
		parser.value_of_key("optimizer").empty() ? 
		"sgd" : parser.value_of_key("optimizer");
	
	const int memory_window = 
		parser.value_of_key("memory_window").empty() ?
		30 : std::atoi(parser.value_of_key("memory_window").c_str());

	const int num_hidden_layers = 
		parser.value_of_key("num_hidden_layers").empty() ? 
		1 : std::atoi(parser.value_of_key("num_hidden_layers").c_str());
	
	const int seed = 
		parser.value_of_key("seed").empty() ?
		0 : std::atoi(parser.value_of_key("seed").c_str());
	
	const int instances = 
		parser.value_of_key("instances").empty() ?
		100000 : std::atoi(parser.value_of_key("instances").c_str());
	
	const std::string target_activation_type = 
		parser.value_of_key("target_activation").empty() ? 
		"tanh" : parser.value_of_key("target_activation");
	
	const double validation_threshold = 
		parser.value_of_key("validation_threshold").empty() ?
		5e-3 : std::atof(parser.value_of_key("validation_threshold").c_str());

	// required inputs
	const std::string chain = parser.value_of_key("chain") + "/";
	const std::string qubits = parser.value_of_key("qubits");

	// derived inputs
	const std::string input_root = data_root + "/input/";
	const std::string output_root = data_root + "/output/";
	const std::string base_input_dir = input_root + chain;
	const std::string base_output_dir = output_root + chain;
	
	const std::string data_path = base_input_dir + qubits + "-qubits.bin";

	Dataset<double> dataset(std::atoi(qubits.c_str()), data_path, batch_size, instances);

	const int input_layer_size = dataset.feature_length();
	const int output_layer_size = dataset.target_length();

	Activation<double> input_activation("linear");
	Activation<double> hidden_activation(hidden_activation_type);
	Activation<double> target_activation(target_activation_type);

	Topology<double> topology(batch_size);
	topology.push_back(input_layer_size, input_activation);

	for (int i = 0; i < num_hidden_layers; ++i) 
		topology.push_back(hidden_layer_size, hidden_activation);
	
	topology.push_back(output_layer_size, target_activation);

	srand(static_cast<unsigned int>(seed));
	Network<double> network(topology);
	network.validate_topology(dataset);

	std::vector<Network<double>> networks;

	for (int i = 0; i < dataset.num_eigenvectors(); ++i) networks.push_back(network);

	Operators<double> operators(std::atoi(qubits.c_str()));
	Optimizer<double> optimizer(optimizer_type, learning_rate, decay_rate, epsilon_conditioner);

	Loss<double> loss(loss_type, operators, lagrange_multiplier, trade_off_parameter, random_domain_bound);
	Model<double> model(networks, loss, optimizer);

	//std::ofstream mse_stream(base_output_dir + "mse-" + loss_type + ".dat", std::ofstream::trunc);
	boost::circular_buffer<double> mse_history(memory_window);

	for (int i = 0; i < memory_window; ++i)
		mse_history.push_back(std::numeric_limits<double>::max());

	double cutout_value, inactive_value;
	std::vector<double> activities(memory_window);

	//std::cout << std::scientific;
	
	//mse_stream << "#Epoch   Training      Testing      Validation\n";
	//mse_stream << std::scientific;

	for (int epoch = 0; epoch <= epochs; ++epoch)
	{
		mse_history.push_back(model.mse(dataset));

		std::adjacent_difference(
				mse_history.begin(), 
				mse_history.end(), 
				activities.data(),
				[](double x, double y){return std::abs(x-y);});

		inactive_value = std::accumulate(activities.begin() + 1, activities.end(), 0.0);
		inactive_value /= memory_window;

		cutout_value = std::accumulate(mse_history.begin(), mse_history.end(), 0.0);
		cutout_value /= memory_window;

		if (cutout_option == "true" && cutout_value < validation_threshold) 
		{
			open_ascii_escape("green");
			std::cout << "\nSUCCESS! MODEL SATISFIED VALIDATION CUTOUT\n" << std::endl;
			close_ascii_escape();
			break;
		}
		else if (cutout_option == "true" && inactive_value < inactive_threshold)
		{
			open_ascii_escape("yellow");
			std::cout << "\nMAJOR ISSUE! MODEL STAGNATED\n" << std::endl;
			close_ascii_escape();
			break;

		}
		else
		{
			//std::cout << epoch << '\t' << training_mse << '\t'
				//<< testing_mse << '\t' << mse_history.back() << std::endl;
			
			//mse_stream<< epoch << '\t' << training_mse << '\t'
				//<< testing_mse << '\t' << mse_history.back() << std::endl;

		
			//if (parser.value_of_key("write_average_magnetization") == "true")
				//model.write_average_magnetization(testing_data, 
						//base_output_dir + "testing/avg-" + loss_type, epoch);

			//if (parser.value_of_key("write_coefficients") == "true")
				//model.write_coefficients(testing_data, 
						//base_output_dir + "testing/coe-" + loss_type, epoch);

			//if (parser.value_of_key("write_lyapunov_estimate") == "true") 
				//model.write_lyapunov_estimate(testing_data, 
						//base_output_dir + "testing/lya-" + loss_type, epoch);

			//if (parser.value_of_key("write_schrodinger_error") == "true")
				//model.write_schrodinger_error(testing_data, 
						//base_output_dir + "testing/sch-" + loss_type, epoch);
						
			//std::cout << "#Epoch      MSE  \n";
			
			std::cout << epoch << '\t' << mse_history.back() << std::endl;


			model.learn_from(dataset);

			//if (parser.value_of_key("shuffle") == "true")
				//training_data.shuffle();
		}

	}

	model.predictive_power(dataset, 0);

	std::cout << std::scientific;
	std::cout << lagrange_multiplier << '\t' << model.pure_cost(dataset) << '\n';

	//if (parser.value_of_key("save_model") == "true") 
		//model.save(base_output_dir + loss_type + "-model.net");


	return 0;
}
