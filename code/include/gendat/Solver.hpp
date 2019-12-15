///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Solver.hpp ***                              //
//                                                                           //
// created December 13, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Solver_hpp
#define Solver_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

template<typename T>
using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
class Solver : public Eigen::SelfAdjointEigenSolver<DenseMatrix<T>>
{
	public:
		Solver(void) : Eigen::SelfAdjointEigenSolver<DenseMatrix<T>>() {}

		Solver& operator=(const Eigen::SelfAdjointEigenSolver<DenseMatrix<T>>& source)
		{
			this->Eigen::SelfAdjointEigenSolver<DenseMatrix<T>>::operator=(source);
			return *this;
		}

		void append_to_file(std::string fpath) const;
		void read(std::string fpath, std::streampos pos = std::ios::beg);
		void print(void) const;
};


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
	
#endif /* Solver_hpp */

