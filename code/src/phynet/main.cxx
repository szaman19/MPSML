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
#include <boost/circular_buffer.hpp>
#include <phynet/Model.hpp>
#include <phynet/Parser.hpp>
#include <phynet/Dataset.hpp>
#include <string>

//int main()
//{
	//int num_qubits = 4;
	//int batch_size = 100;

	//std::string N = std::to_string(num_qubits);

	//std::string root = "/home/csingh5/Documents/guided-machine-learning/code/data/";
	//std::string fpath = root + "input/clean-Ising/" + N + "-qubits.bin";

	//Dataset<double> data(num_qubits, fpath, batch_size);

	//return 0;
//}

int main( int argc, char *argv[] )
{
	std::string param_file;

	if (argc == 2) param_file = argv[1];
	else param_file = "etc/phynet.in";

	Parser parser(param_file);
	
	const std::string loss_type = 
		parser.value_of_key("loss_type").empty() ? 
		"bb" : parser.value_of_key("loss_type");
	
	const std::size_t batch_size = 
		parser.value_of_key("batch_size").empty() ? 
		100 : std::atoi(parser.value_of_key("batch_size").c_str());

	const std::string cutout_option = 
		parser.value_of_key("cutout_option").empty() ? 
		"true" : parser.value_of_key("cutout_option");
	
	const std::size_t epochs = 
		parser.value_of_key("epochs").empty() ? 
		1000 : std::atoi(parser.value_of_key("epochs").c_str());
	
	const std::string data_root = 
		parser.value_of_key("data_root").empty() ? 
		"data" : parser.value_of_key("data_root");

	const std::string hidden_activation_type = 
		parser.value_of_key("hidden_activation").empty() ? 
		"tanh" : parser.value_of_key("hidden_activation");
	
	const std::size_t hidden_layer_size =
		parser.value_of_key("hidden_layer_size").empty() ? 
		100 : std::atoi(parser.value_of_key("hidden_layer_size").c_str());
	
	const double inactive_threshold = 
		parser.value_of_key("inactive_threshold").empty() ?
		1e-6 : std::atof(parser.value_of_key("inactive_threshold").c_str());

	const double lagrange_multiplier = 
		parser.value_of_key("lagrange_multiplier").empty() ?
		0.1 : std::atof(parser.value_of_key("lagrange_multiplier").c_str());
	
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
	
	const std::size_t memory_window = 
		parser.value_of_key("memory_window").empty() ?
		30 : std::atoi(parser.value_of_key("memory_window").c_str());

	const std::size_t num_hidden_layers = 
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

	std::cout << data_path << '\n';

	Dataset<float> dataset(std::atoi(qubits.c_str()), data_path, batch_size, instances);

	// technically shouldn't have to run this everytime, but its harmless and guarantees existance 
	//training_data.generate_template_average_file(base_output_dir+"training/avg-sz-exact.dat");
	//testing_data.generate_template_average_file(base_output_dir+"testing/avg-sz-exact.dat");
	//validation_data.generate_template_average_file(base_output_dir+"validation/avg-sz-exact.dat");
	
	const std::size_t input_layer_size = dataset.feature_length();
	const std::size_t output_layer_size = dataset.target_length();

	Activation<float> input_activation("linear");
	Activation<float> hidden_activation(hidden_activation_type);
	Activation<float> target_activation(target_activation_type);

	Topology<float> topology(batch_size);
	topology.push_back(input_layer_size, input_activation);

	for (std::size_t i = 0; i < num_hidden_layers; ++i) 
		topology.push_back(hidden_layer_size, hidden_activation);
	
	topology.push_back(output_layer_size, target_activation);

	srand(static_cast<unsigned int>(seed));
	Network<float> network(topology);
	network.validate_topology(dataset);

	std::vector<Network<float>> networks;

	for (int i = 0; i < dataset.num_eigenvectors(); ++i) networks.push_back(network);

	Loss<float> loss(loss_type, lagrange_multiplier, trade_off_parameter);
	Optimizer<float> optimizer(optimizer_type, learning_rate, decay_rate, epsilon_conditioner);

	Model<float> model(networks, loss, optimizer);

	//std::ofstream mse_stream(base_output_dir + "mse-" + loss_type + ".dat", std::ofstream::trunc);
	//boost::circular_buffer<double> mse_history(memory_window);

	//for (std::size_t i = 0; i < memory_window; ++i)
		//mse_history.push_back(std::numeric_limits<double>::max());

	//double cutout_value, training_mse, testing_mse, inactive_value;
	//std::vector<double> activities(memory_window);

	//std::cout << "#Epoch   Training      Testing      Validation\n";
	//std::cout << std::scientific;
	
	//mse_stream << "#Epoch   Training      Testing      Validation\n";
	//mse_stream << std::scientific;

	//std::cout << "Epoch \t MSE \t PredP \n";

	double cutout_value = 1e6;
	double inactive_value = 1e6;

	for (std::size_t epoch = 0; epoch <= epochs; ++epoch)
	{
		//mse_history.push_back(model.mse(validation_data));
		//training_mse = model.mse(training_data);
		//testing_mse = model.mse(testing_data);

		//std::adjacent_difference(
				//mse_history.begin(), 
				//mse_history.end(), 
				//activities.data(),
				//[](double x, double y){return std::abs(x-y);});

		//inactive_value = std::accumulate(activities.begin() + 1, activities.end(), 0.0);
		//inactive_value /= memory_window;

		//cutout_value = std::accumulate(mse_history.begin(), mse_history.end(), 0.0);
		//cutout_value /= memory_window;

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
						
			std::cout << epoch << '\t' << model.mse(dataset) 
					  << '\t' << model.predictive_power(dataset, epoch) << std::endl;

			model.learn_from(dataset);

			//if (parser.value_of_key("shuffle") == "true")
				//training_data.shuffle();
		}

	}

	//if (parser.value_of_key("save_model") == "true") 
		//model.save(base_output_dir + loss_type + "-model.net");


	return 0;
}

	//af::setDevice(1);
	//af::setBackend(AF_BACKEND_CPU);
	//af::info();

	//const int N = 1000;
	//const int trials = 1;

	//eigen e1(N,N), e2(N,N);
	//e1.setRandom();

	//af::timer::start();
	//for (int i = 0; i < trials; ++i) e2 = e1 * e1;
	//std::cout << "CPU runtime : " << af::timer::stop()/trials << " seconds\n";

	//// Select a device and display arrayfire info
	//af::info();

	//af::array a1(N, N, e1.data()), a2(N,N);
	//af::timer::start();
	//for (int i = 0; i < trials; ++i) a2 = matmul(a1,a1);
	//std::cout << "GPU runtime : " << af::timer::stop()/trials << " seconds\n";
	   
	//Eigen::Map<Eigen::MatrixXf> res(a2.host<double>(), N, N);
	//if (res.isApprox(e2)) std::cout << "CPU and GPU agree\n";
	//else std::cout << "CPU and GPU disagree\n";
