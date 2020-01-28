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
Generator<T>::Generator(std::string model, 
		int num_qubits, int num_transverse, int num_realizations)
	: 
	num_qubits(num_qubits), 
	num_transverse(num_transverse),
	num_realizations(num_realizations),
	model(model)
{
	if (model != "ising" && model != "xyz")
	{
		std::cerr << "UNSUPPORTED MODEL" << std::endl;
		exit(-1);
	}
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
	Eigen::Array<T, Eigen::Dynamic, 1> Bx;

	if (model == "ising")
		Bx = Eigen::Array<T, Eigen::Dynamic, 1>::LinSpaced(num_transverse, BX_MIN, BX_MAX);

	if (model == "xyz")
	{
		std::vector<T> x1 {-0.61, -0.65, -0.67, -0.56, -0.43, -0.22, -0.08, 
			0.12, 0.45, 0.86, 1.14, 1.35, 1.45, 1.52, 1.50};

		std::vector<T> y1 {-0.54, -0.35, -0.16, -0.11,  0.32,  0.55,  0.67, 
			0.82, 1.07, 1.15, 1.06, 0.85, 0.67, 0.42, 0.23};


		Bx = Eigen::Array<T, Eigen::Dynamic, 1>::LinSpaced(num_transverse, -2, 2);
	}

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

			if (model == "ising")
				solver.compute(operators.ising_hamiltonian(fields));

			if (model == "xyz")
				solver.compute(operators.xyz_hamiltonian(fields));

			if (solver.info() != Eigen::Success)
			{
				std::cerr << " \nDIAGONALIZATION FAILED\n\n";
				exit(-1);
			}

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
