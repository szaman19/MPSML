#define EIGEN_USE_MKL_ALL
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/filesystem.hpp>
#include <vector>
#include <boost/program_options.hpp>

using namespace std;
using namespace Eigen;

#define NUM_PARAMS_PER_LINE 5

int main(int argc, char *argv[])
{
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	desc.add_options()
		("help", "display help message")
		("fpath", po::value<std::string>(), "path to stats file")
		("algo", po::value<std::string>(), "choost bb, pv, mpv")
		("qubits", po::value<std::string>(), "set number of qubits");
		
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		exit(-1);
	}

	if (!vm.count("fpath") || !vm.count("algo") || !vm.count("qubits")) 
	{
		std::cout << "Must specify at least fpath, algo, qubits\n";
		std::cout << "Please run with --help to set relavant parameters" << std::endl;
		exit(-1);
	}

	std::string fpath = vm["fpath"].as<std::string>();
	std::string algo = vm["algo"].as<std::string>();
	std::string qubits = vm["qubits"].as<std::string>();

	long bytes_per_line = NUM_PARAMS_PER_LINE*sizeof(double);
	long lines = boost::filesystem::file_size(fpath)/bytes_per_line;

	ifstream file(fpath, std::ios::in | std::ios::binary);

	Array<double, Dynamic, Dynamic, RowMajor> data(lines, NUM_PARAMS_PER_LINE);
	file.read((char*)data.data(), data.size()*sizeof(double));

	int trials = data.col(0).maxCoeff() + 1;
	int epochs = data.col(1).maxCoeff() + 1;
	int instances = lines / (trials * epochs);

	//cout << "Lines:   " << lines << '\n';
	//cout << "Trials:  " << trials << '\n';
	//cout << "Epochs:  " << epochs << '\n';
	//cout << "Samples: " << instances << '\n';

	Array<double, Dynamic, Dynamic, RowMajor> average_as_function_of_epochs(epochs, 2);
	Array<double, Dynamic, Dynamic, RowMajor> variance_as_function_of_epochs(epochs, 2);
	Array<double, Dynamic, Dynamic, RowMajor> lyapunov_as_function_of_bx(instances, 1);

	average_as_function_of_epochs.setZero();
	variance_as_function_of_epochs.setZero();

	for (int line = 0; line < lines; ++line)
	{
		average_as_function_of_epochs((int)data(line, 1), 0) = data(line, 1) + 1;
		average_as_function_of_epochs((int)data(line, 1), 1) += fabs(data(line, 4)) / (lines/epochs);
	}
	
	for (int line = 0; line < lines; ++line)
	{
		variance_as_function_of_epochs((int)data(line, 1), 0) += data(line, 4)*data(line, 4);
		variance_as_function_of_epochs((int)data(line, 1), 1) += fabs(data(line, 4));
	}

	for (int i = 0; i < variance_as_function_of_epochs.rows(); ++i)
	{
		variance_as_function_of_epochs(i, 1) *= variance_as_function_of_epochs(i, 1); 		
		variance_as_function_of_epochs(i, 1) /= (lines/epochs);
	}
	
	for (int i = 0; i < variance_as_function_of_epochs.rows(); ++i)
	{
		variance_as_function_of_epochs(i, 1) = variance_as_function_of_epochs(i, 0) - 
			variance_as_function_of_epochs(i, 1);

		variance_as_function_of_epochs(i, 1) /= ( (lines/epochs) - 1 );
		variance_as_function_of_epochs(i, 0) = i; 		
	}

	std::ofstream ofile(algo + "-" + qubits + "-mean-and-variance-vs-epochs.dat");
	
	ofile << std::scientific;
	for (int i = 0; i < variance_as_function_of_epochs.rows(); ++i)
	{
		ofile << i << '\t';
		ofile << average_as_function_of_epochs(i,1) << '\t';
		ofile << variance_as_function_of_epochs(i,1) << '\n';
	}


	return 0;
}









