///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Generator.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <gendat/Generator.hpp>

template <typename T>
Generator<T>::Generator(int num_qubits, int num_transverse, int num_realizations)
	: 
	num_qubits(num_qubits), 
	num_transverse(num_transverse),
	num_realizations(num_realizations)
{
	// null body
}
template <typename T>
void Generator<T>::set_dump_location(std::string fpath)
{
	//prompt_if_file(fpath);
	this->fpath = fpath;
	if (boost::filesystem::exists(fpath)) boost::filesystem::remove(fpath);
	ready_to_dump = true;
}

template <typename T>
void Generator<T>::run(void) const
{
	Eigen::Array<T, Eigen::Dynamic, 1> Bx = 
		Eigen::Array<T, Eigen::Dynamic, 1>::LinSpaced(num_transverse, BX_MIN, BX_MAX);

	T W;

	if (num_realizations > 1)
		W = DISORDER_STRENGTH;
	else
		W = 0;

	Fields<T> fields;
	Solver<T> solver;
	Operators<T> operators(num_qubits);

	for (int iBx = 0; iBx < Bx.size(); ++iBx)
	{
		for (int iW = 0; iW < num_realizations; ++iW)
		{
			fields = Fields<T>(num_qubits, J, Bx[iBx], BZ, W);
			solver.compute(operators.hamiltonian(fields));

			if (solver.info() != Eigen::Success)
			{
				std::cerr << " \nDIAGONALIZATION FAILED\n\n";
				exit(-1);
			}

			if (solver.eigenvectors().col(0).mean() < 0)
				Instance<T>(fields, solver).append_to_file(fpath);
		}
	}
}

template <typename T>
void Generator<T>::prompt_if_file(std::string fpath) const 
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

template class Generator<double>;
