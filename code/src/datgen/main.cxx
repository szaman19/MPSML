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
#include "operators.hpp"
#include "Decomp.hpp"

int main( /* int argc, char *argv[] */ )
{
	const int num_qubits = 4;
	const int num_disorder_realizations = 10;
	const int num_bx_values = 10;
	
	double W;

	Eigen::ArrayXf bx_values = Eigen::ArrayXf::LinSpaced(num_bx_values, 0.1, 2);

	saes_t s;

	if (num_disorder_realizations > 1) W = 0.01;
	else W = 0.0;

	int drs, bxs;
	drs = 0;
	bxs = 0;

	for (int i = 0; i < bx_values.size(); ++i)
	{
		for (int j = 0; j < num_disorder_realizations; ++j)
		{
			std::cout << "bx: " << bx_values[i] << '\t';
			s.compute(hamiltonian(num_qubits, Fields(num_qubits, bx_values[i], W)));
			if (s.info() != Eigen::Success)
			{
				std::cerr << "failed diag\n";
				exit(-1);
			}

			drs++;
		}

		std::cout << "ran through " << drs << " disorder realizations\n";
		drs = 0;
		bxs++;
	}

	std::cout << "and " << bxs << " Bx values" << '\n';


	return 0;
}


