///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** datgen main.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Sparse>
#include <gendat/Fields.hpp>
#include <gendat/Generator.hpp>
#include <gendat/Operators.hpp>
#include <boost/program_options.hpp>

int main(int argc, char *argv[])
{
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	desc.add_options()
		("help", "display help message")
		("model", po::value<std::string>(), "choose ising or xyz")
		("qubits", po::value<int>(), "set number of qubits")
		("dBx", po::value<int>(), "partition along dimension Bx")
		("dBz", po::value<int>(), "partition along dimension Bz")
		("Bx_min", po::value<double>(), "min transverse field")
		("Bx_max", po::value<double>(), "max transverse field")
		("Bz_min", po::value<double>(), "min longitudinal field")
		("Bz_max", po::value<double>(), "max longitudinal field")
		("fpath", po::value<std::string>(), "location to write data")
		("replicas", po::value<int>(), "set number of disorder replicas")
		("disorder", po::value<double>(), "set disorder strength")
		("coupling", po::value<double>(), "Ising coupling J");
		

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		exit(-1);
	}

	if (!vm.count("qubits") || !vm.count("fpath")) 
	{
		std::cout << "Must specify at least qubits and fpath\n";
		std::cout << "Please run with --help to set relavant parameters" << std::endl;
		exit(-1);
	}

	std::string model = vm.count("model") ? vm["model"].as<std::string>() : "ising";

	Generator<double>generator(model);	

	// required inputs 
	generator.qubits = vm["qubits"].as<int>();
	generator.fpath = vm["fpath"].as<std::string>();

	generator.replicas = vm.count("replicas") ? vm["replicas"].as<int>() : 1;
	generator.coupling = vm.count("coupling") ? vm["coupling"].as<double>() : 1;
	generator.disorder = vm.count("disorder") ? vm["disorder"].as<double>() : 0.01;
	
	generator.dBx = vm.count("dBx") ? vm["dBx"].as<int>() : 10;
	generator.dBz = vm.count("dBz") ? vm["dBz"].as<int>() : 10;

	generator.Bx_min = vm.count("Bx_min") ? vm["Bx_min"].as<double>() : 0.0;
	generator.Bx_max = vm.count("Bx_max") ? vm["Bx_max"].as<double>() : 2.0;
	
	generator.Bz_min = vm.count("Bz_min") ? vm["Bz_min"].as<double>() : 0.0;
	generator.Bz_max = vm.count("Bz_max") ? vm["Bz_max"].as<double>() : 2.0;

	generator.run();

	return 0;
}


