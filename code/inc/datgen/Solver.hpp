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

class Solver : public Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
{
	public:
		Solver(void) : Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>() {}

		Solver& operator=(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& source)
		{
			this->Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>::operator=(source);
			return *this;
		}

		void append_to_file(std::string fpath) const;
		void read(std::string fpath, std::streampos pos = std::ios::beg);
		void print(void) const;

		void reserve(int num_qubits)
		{
			int dim = (int)(pow(2, num_qubits) + 0.5);
			this->compute(Eigen::MatrixXd::Identity(dim, dim));
		}
};


	
#endif /* Solver_hpp */

