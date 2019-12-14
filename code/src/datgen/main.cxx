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
#include "Decomp.hpp"
#include "Operators.hpp"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		if (argv[0])
		{
			std::cout << " Usage: " << argv[0];
			std::cout << " <num_qubits> ";
			std::cout << " <num_transverse_fields> ";
			std::cout << " <num_disorder_realizations> ";
			std::cout << " <fpath>";
			std::cout << std::endl;
		}
		else
		{
			std::cout << " Usage: datgen.x";
			std::cout << " <num_qubits> ";
			std::cout << " <num_transverse_fields> ";
			std::cout << " <num_disorder_realizations> ";
			std::cout << " <fpath>";
			std::cout << std::endl;
		}

		exit(-1);
	}

	std::stringstream a1(argv[1]);
	std::stringstream a2(argv[2]);
	std::stringstream a3(argv[3]);

	int num_qubits, num_transverse_fields, num_disorder_realizations;

	a1 >> num_qubits;
	a2 >> num_transverse_fields;
	a3 >> num_disorder_realizations;

	Decomp decomp(num_qubits, num_transverse_fields, num_disorder_realizations);	

	decomp.compute();

	decomp.write(argv[4]);

	return 0;
}


