///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Generator.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Generator.hpp>

template <typename T> 
Generator<T>::Generator(){}

template <typename T>
void Generator<T>::run(void) const
{
	// Remove boost dependecies 
	// if (boost::filesystem::exists(fpath)) boost::filesystem::remove(fpath); 

	Eigen::Array<T, Eigen::Dynamic, 1> Bx, Bz;

	Bx = Eigen::Array<T, Eigen::Dynamic, 1>::LinSpaced(dBx, Bx_min, Bx_max);
	Bz = Eigen::Array<T, Eigen::Dynamic, 1>::LinSpaced(dBz, Bz_min, Bz_max);

	if (dBz == 1) Bz(0) = 0.01; // Why?
	if (dBx == 1) Bx(0) = 0.01; // Why?

	T W = 0;

	Fields<T> fields;
	Solver<T> solver;
	Operators<T> operators(qubits);
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop, origin;
	std::chrono::duration<double> elapsed;
	origin = std::chrono::high_resolution_clock::now();
	elapsed = origin - origin;

	for (int iBx = 0; iBx < Bx.size(); ++iBx)
	{
		for (int iBz = 0; iBz < Bz.size(); ++iBz)
		{
			for (int iW = 0; iW < 1; ++iW)
			{
				fields = Fields<T>(qubits, coupling, Bx[iBx], Bz[iBz], W);

				start = std::chrono::high_resolution_clock::now();
				solver.compute(operators.ising_hamiltonian(fields));
				stop = std::chrono::high_resolution_clock::now();
				elapsed += stop - start;

				if (solver.info() != Eigen::Success)
				{
					std::cerr << " \nDIAGONALIZATION FAILED\n\n";
					exit(-1);
				}

				Instance<T>(fields, solver).append_to_file(fpath);
			}
		}
	}
	std::cout << "Diagonalization time: " << elapsed.count() << " seconds\n";
}

template class Generator<double>;
