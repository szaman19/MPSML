///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Solver.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Solver.hpp"

void Solver::append_to_file(std::string fpath) const
{
	std::ofstream file(fpath.c_str(), std::ios::app | std::ios::binary);

	file.write((char*)this->eigenvalues().data(), this->eigenvalues().size() * sizeof(double));
	file.write((char*)this->eigenvectors().data(), this->eigenvectors().size() * sizeof(double));
}

void Solver::read(std::string fpath, std::streampos pos)
{
	std::ifstream file(fpath.c_str(), std::ios::binary);

	if (file.is_open())
	{
		file.seekg(pos);
		file.read((char*)this->eigenvalues().data(), this->eigenvalues().size()*sizeof(double));
		file.read((char*)this->eigenvectors().data(), this->eigenvectors().size()*sizeof(double));
	}
	else
	{
		std::cerr << "ERROR: FAILED TO OPEN " << fpath << '\n';
		exit(-1);
	}
}

void Solver::print(void) const 
{
	std::cout << "Eigenvalues:\n";
	std::cout << this->eigenvalues().transpose() << "\n\n";
	std::cout << "Eigenvectors:\n";
	std::cout << this->eigenvectors() << "\n\n";
}
