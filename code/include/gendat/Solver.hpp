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
#include <math.h>

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


	
#endif /* Solver_hpp */

