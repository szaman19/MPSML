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

int main(int argc, char *argv[])
{
	if (argc != 6)
	{
		if (argv[0])
		{
			std::cout << " Usage: " << argv[0];
			std::cout << " <ising or xyz> ";
			std::cout << " <num_qubits> ";
			std::cout << " <num_parameter_points> ";
			std::cout << " <num_disorder_realizations> ";
			std::cout << " <fpath>";
			std::cout << std::endl;
		}
		else
		{
			std::cout << " Usage: gendat";
			std::cout << " <ising or xyz> ";
			std::cout << " <num_qubits> ";
			std::cout << " <num_parameter_points> ";
			std::cout << " <num_disorder_realizations> ";
			std::cout << " <fpath>";
			std::cout << std::endl;
		}

		exit(-1);
	}

	std::stringstream a1(argv[1]);
	std::stringstream a2(argv[2]);
	std::stringstream a3(argv[3]);
	std::stringstream a4(argv[4]);

	int num_qubits, num_transverse_fields, num_disorder_realizations;
	std::string model;

	a1 >> model;
	a2 >> num_qubits;
	a3 >> num_transverse_fields;
	a4 >> num_disorder_realizations;

	Generator<double> 
		generator(model, num_qubits, num_transverse_fields, num_disorder_realizations);	

	generator.set_dump_location(argv[5]);

	generator.run();

	return 0;
}


