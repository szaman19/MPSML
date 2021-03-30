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
#include <chrono>
#include <boost/circular_buffer.hpp>
#include <boost/progress.hpp>
#include <phynet/Model.hpp>
#include <phynet/Parser.hpp>
#include <phynet/Dataset.hpp>
#include <gendat/Operators.hpp>
#include <phynet/Optimizer.hpp>

int main( int argc, char *argv[] )
{
	std::cout << std::endl;
	std::string param_file;

	if (argc == 2) param_file = argv[1];
	else param_file = "etc/phynet.in";

	Parser parser(param_file);
	
	// #################### GATHER USER INPUTS ################### // 
	// required
	const std::string chain = parser.value_of_key("chain") + "/";
	const std::string phase = parser.value_of_key("phase");
	const std::string qubits = parser.value_of_key("qubits");

	// defaults provided
	const std::string loss_type = 
		parser.value_of_key("loss_type").empty() ? 
		"bb" : parser.value_of_key("loss_type");
	
	const int batch_size = 
		parser.value_of_key("batch_size").empty() ? 
		10 : std::atoi(parser.value_of_key("batch_size").c_str());

	const std::string cutout_option = 
		parser.value_of_key("cutout_option").empty() ? 
		"true" : parser.value_of_key("cutout_option");
	
	const int epochs = 
		parser.value_of_key("epochs").empty() ? 
		10000 : std::atoi(parser.value_of_key("epochs").c_str());
	
	const std::string data_root = 
		parser.value_of_key("data_root").empty() ? 
		"data" : parser.value_of_key("data_root");
	
	const std::string input = 
		parser.value_of_key("input").empty() ? 
		"fields" : parser.value_of_key("input");

	const std::string hidden_activation_type = 
		parser.value_of_key("hidden_activation").empty() ? 
		"tanh" : parser.value_of_key("hidden_activation");
	
	const std::string hidden_layer_dimensions = 
		parser.value_of_key("hidden_layer_dimensions").empty() ? 
		" 10 " : parser.value_of_key("hidden_layer_dimensions");
	
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
		0.01 : std::atof(parser.value_of_key("learning_rate").c_str());
	
	const std::string optimizer_type = 
		parser.value_of_key("optimizer").empty() ? 
		"sgd" : parser.value_of_key("optimizer");
	
	const int memory_window = 
		parser.value_of_key("memory_window").empty() ?
		30 : std::atoi(parser.value_of_key("memory_window").c_str());
	
	const int num_eigenvectors = 
		parser.value_of_key("num_eigenvectors").empty() ?
		1 : std::atoi(parser.value_of_key("num_eigenvectors").c_str());

	const int seed = 
		parser.value_of_key("seed").empty() ?
		rand() : std::atoi(parser.value_of_key("seed").c_str());
	
	const int statistical_trials = 
		parser.value_of_key("statistical_trials").empty() ?
		1 : std::atoi(parser.value_of_key("statistical_trials").c_str());
	
	const std::string target_activation_type = 
		parser.value_of_key("target_activation").empty() ? 
		"tanh" : parser.value_of_key("target_activation");
	
	const std::string post_process = 
		parser.value_of_key("post_process").empty() ? 
		"false" : parser.value_of_key("post_process");

	const double validation_threshold = 
		parser.value_of_key("validation_threshold").empty() ?
		1e-4 : std::atof(parser.value_of_key("validation_threshold").c_str());

	const double epochs_per_cycle = 
		parser.value_of_key("validation_threshold").empty() ?
		epochs/10 : std::atof(parser.value_of_key("validation_threshold").c_str());

	const std::string verbose = 
		parser.value_of_key("verbose").empty() ? 
		"false" : parser.value_of_key("verbose").c_str();
	
	const std::string learning_rate_schedule = 
		parser.value_of_key("learning_rate_schedule").empty() ? 
		"constant" : parser.value_of_key("learning_rate_schedule").c_str();

	// #################### CREATE DATASET ################### // 
	const std::string input_root = data_root + "/input/";
	const std::string output_root = data_root + "/output/";
	const std::string base_input_dir = input_root + chain;
	const std::string base_output_dir = output_root + chain;
	const std::string stats_fpath = base_output_dir + qubits + "-" + loss_type + "-stats.bin";

	const std::string data_path = base_input_dir + qubits + "-qubits.bin";

	Operators<double> operators(std::atoi(qubits.c_str()));
	Dataset<double> dataset(std::atoi(qubits.c_str()), data_path, 
			batch_size, input, operators, phase);

	const int input_layer_size = dataset.feature_length();
	const int output_layer_size = dataset.target_length();

	Activation<double> input_activation("linear");
	Activation<double> hidden_activation(hidden_activation_type);
	Activation<double> target_activation(target_activation_type);

	// #################### CREATE TOPOLOGY ################### // 
	Topology<double> topology(batch_size);

	// set first layer
	topology.push_back(input_layer_size, input_activation);

	std::istringstream layer_dims(hidden_layer_dimensions);

	int tmp;
	std::vector<int> v;
	while(layer_dims >> tmp) v.push_back(tmp);

	// set hidden layers 
	for (auto i : v) topology.push_back(i, hidden_activation);
	
	// set last layer 
	topology.push_back(output_layer_size, target_activation);
	
	std::cout << "Layers: " << v.size() + 2 << '\n';
	std::cout << "Feature length: " << input_layer_size << '\n';
	std::cout << "Target length: " << output_layer_size << '\n';
	std::cout << "Hidden Topology: ";
	for (auto i : v) std::cout << i << '\t';
	std::cout << std::endl;
	
	// #################### INSTANTIATE NETWORKS ################### // 
	srand(static_cast<unsigned int>(seed));
	Network<double> network(topology);
	network.validate_topology(dataset);

	std::vector<Network<double>> networks;
	for (int i = 0; i < num_eigenvectors; ++i) networks.push_back(network);

	Optimizer<double> optimizer(optimizer_type, learning_rate_schedule);
	optimizer.set_learning_rate_min(learning_rate);
	optimizer.set_learning_rate_max(learning_rate*10);
	optimizer.set_epochs_per_cycle(epochs_per_cycle);
	optimizer.set_decay_rate(decay_rate);
	optimizer.set_epsilon_conditioner(epsilon_conditioner);

	Loss<double> loss(loss_type, operators, lagrange_multiplier, 
			trade_off_parameter, random_domain_bound);

	Model<double> model(networks, loss, optimizer);

	// #################### MAIN TRAINING LOOP ################### // 
	boost::circular_buffer<double> mse_history(memory_window);

	for (int i = 0; i < memory_window; ++i)
		mse_history.push_back(std::numeric_limits<double>::max());

	double cutout_value, inactive_value;
	std::vector<double> activities(memory_window);

	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop, origin;
	std::chrono::duration<double> elapsed;
	origin = std::chrono::high_resolution_clock::now();
	elapsed = origin - origin;

	for (int trial = 0; trial < statistical_trials; ++trial)
	{
		std::cout << "\n############ Trial " << trial << " ############\n";
		std::cout << "\nEpochs Exhausted:";
		boost::progress_display progress(epochs);
		for (int epoch = 0; epoch < epochs; ++epoch)
		{
			if (post_process == "true") model.append_wandb_for_radviz("wandb-radviz.bin");

			if (statistical_trials > 1)
				model.append_metrics(dataset, trial, epoch, stats_fpath);

			optimizer.set_epoch(epoch);
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
				std::cout << "Epochs used: " << epoch << '\n';
				break;
			}
			else if (cutout_option == "true" && inactive_value < inactive_threshold)
			{
				open_ascii_escape("yellow");
				std::cout << "\nISSUE! MODEL STAGNATED\n" << std::endl;
				close_ascii_escape();
				std::cout << "Epochs used: " << epoch << '\n';
				break;
			}
			else
			{
				if (verbose == "true") 
					std::cout << epoch << '\t' << mse_history.back() << std::endl;

				start = std::chrono::high_resolution_clock::now();
				model.learn_from(dataset);
				stop = std::chrono::high_resolution_clock::now();
				elapsed += stop - start;
			}

			open_ascii_escape("cyan");
			++progress;
			close_ascii_escape();

			if (cutout_option == "true" && epoch == epochs-1) 
			{
				open_ascii_escape("red");
				std::cout << "WARNING! EXHAUSTED ALL EPOCHS\n" << std::endl;
				close_ascii_escape();
				std::cout << "Epochs used: " << epoch+1 << '\n';
			}
		}
		std::cout << "Training time: " << elapsed.count() << " seconds\n";
		std::cout << "Final MSE: " << mse_history.back() << '\n';

		model.print_inference_time(dataset); 
		model.print_average_overlap(dataset);
		model.print_average_sz_error(dataset);

		if (post_process == "true" && num_eigenvectors == 1)
		{
			model.write_magnetization(dataset, "mag.dat");
			model.write_overlap(dataset, "ovr.dat");
			model.write_radial_visualization(dataset, "rad.dat");
			model.write_entanglement_entropy(dataset, "ent.dat");

			std::system("plot-metrics.gnu");
			boost::filesystem::remove("mag.dat");
			boost::filesystem::remove("ent.dat");
			boost::filesystem::remove("ovr.dat");
			boost::filesystem::remove("rad.dat");
			boost::filesystem::remove("anchors.dat");
		}

		model.reset();
		dataset.import();
	}

	std::string ccom = "calc-stats-vs-epochs --fpath " + stats_fpath;
	ccom += " --algo " + loss_type + " --qubits " + qubits;

	std::system(ccom.c_str());

	std::cout << std::endl;
	return 0;
}
