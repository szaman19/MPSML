///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Solver.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Solver.hpp>

template<typename T>
void Solver<T>::append_to_file(std::string fpath) const
{
	std::ofstream file(fpath.c_str(), std::ios::app | std::ios::binary);

	DenseMatrix<T> vecs = this->eigenvectors();

	if (vecs(0,0) < 0 ) vecs *= -1;

	file.write((char*)this->eigenvalues().data(), this->eigenvalues().size() * sizeof(T));
	file.write((char*)vecs.data(), vecs.size() * sizeof(T));
}

template<typename T>
void Solver<T>::print(void) const 
{
	std::cout << "Eigenvalues:\n";
	std::cout << this->eigenvalues().transpose() << "\n\n";
	std::cout << "Eigenvectors:\n";
	std::cout << this->eigenvectors() << "\n\n";
}


template class Solver<double>;
