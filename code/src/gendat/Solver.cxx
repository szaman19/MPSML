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

	DenseMatrix<T> vecs = this->eigenvectors();

	if (vecs(0,0) < 0 ) vecs *= -1;
	
	double eig_val[1] = {this->eigenvalues().data()[0]};
	file.write((char*)eig_val, sizeof(T));
	file.write((char*)vecs.row(0).data(), vecs.row(0).size() * sizeof(T));
	

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
