///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** Decomp.cxx ***                             //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Decomp.hpp"

Decomp::Decomp(int num_qubits, int num_transverse, int num_realizations)
	: 
	num_qubits(num_qubits), 
	num_transverse(num_transverse),
	num_realizations(num_realizations)
{
	// null body
}

void Decomp::compute(void)
{
	Eigen::ArrayXf Bx = 
		Eigen::ArrayXf::LinSpaced(num_transverse, BX_MIN, BX_MAX);

	double W;

	if (num_realizations > 1) 
		W = 0.01;
	else
		W = 0.0;

	Solver waves;
	Fields field;
	Operators operators(num_qubits);

	for (int iBx = 0; iBx < Bx.size(); ++iBx)
	{
		for (int iW = 0; iW < num_realizations; ++iW)
		{
			field = Fields(num_qubits, J, Bx[iBx], BZ, W);
			waves.compute(operators.hamiltonian(field));
			
			if (waves.info() != Eigen::Success)
			{
				std::cerr << " \nDIAGONALIZATION FAILED\n\n";
				exit(-1);
			}
			
			instances.push_back(Instance(field, waves));
		}
	}
}

const Eigen::MatrixXd& Decomp::unitary(int instance) const 
{
	return instances[instance].unitary();
}

const Eigen::VectorXd& Decomp::spectrum(int instance) const 
{
	return instances[instance].spectrum();
}

void Decomp::prompt_if_file(std::string fpath) const 
{
	if (boost::filesystem::exists(fpath))
	{
		std::cout << " WARNING: FILE EXISTS... REMOVE " 
				  << fpath << " [y/n] \n";

		char choice;
		std::cin >> choice;

		if (choice == 'y')
		{
			std::cout << " OK, REMOVING " << fpath << '\n';
			boost::filesystem::remove(fpath);
		}
		else if (choice == 'n') 
		{
			std::cout << " QUIT OR APPEND ANYWAY? [q/a]\n";
			std::cin >> choice;

			if (choice == 'q')
			{
				std::cout << " OK, QUITTING\n";
				exit(-1);
			}
			else if (choice == 'a')
			{
				std::cout << "CONTINUING\n";
			}
		}
	}
}

void Decomp::write(std::string fpath) const 
{
	if (boost::filesystem::exists(fpath)) 
		boost::filesystem::remove(fpath);

	for (std::size_t i = 0; i < instances.size(); ++i)
		instances[i].append_to_file(fpath);	
}

void Decomp::read(std::string fpath)
{
	int dim = (int)(pow(2, num_qubits) + 0.5);
	int off = (4 + num_qubits + dim + dim*dim)*sizeof(double);

	std::streampos pos;
	for (int i = 0; i < num_transverse*num_realizations; ++i)
	{
		pos = off * i;
		instances[i].read(fpath, pos);
	}
}

void Decomp::print(int instance) const
{
	instances[instance].print();
}
