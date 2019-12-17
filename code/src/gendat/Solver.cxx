///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Solver.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <gendat/Solver.hpp>

template<typename T>
void Solver<T>::append_to_file(std::string fpath) const
{
	std::ofstream file(fpath.c_str(), std::ios::app | std::ios::binary);

	file.write((char*)this->eigenvalues().data(), this->eigenvalues().size() * sizeof(T));
	file.write((char*)this->eigenvectors().data(), this->eigenvectors().size() * sizeof(T));
}

template<typename T>
void Solver<T>::print(void) const 
{
	std::cout << "Eigenvalues:\n";
	std::cout << this->eigenvalues().transpose() << "\n\n";
	std::cout << "Eigenvectors:\n";
	std::cout << this->eigenvectors() << "\n\n";
}


template class Solver<float>;
template class Solver<double>;
