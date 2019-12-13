///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** operators.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef operators_hpp
#define operators_hpp
#include <iostream> 
#include <vector>
#include <Eigen/KroneckerProduct>
#include <random>
#include "pauli.hpp"

class Fields
{
	public: 
		Fields(int num_qubits, double Bx, double W = 0, double Bz = 0.05, double J = 1.0)
			: transverse(Bx), disorder_strength(W), longitudinal(Bz), coupling(J)
		{
			std::random_device rd;
			std::default_random_engine engine(rd());
			std::uniform_real_distribution<double> dist(-W, W);

			for (int i = 0; i < num_qubits; ++i) local_fields.push_back(dist(engine));
			for (int i = 0; i < num_qubits; ++i) std::cout << local_fields[i] << '\t';
			std::cout << '\n';
		}

		double transverse;
		double disorder_strength;
		double longitudinal;
		double coupling;

		std::vector<double> local_fields;
};

Eigen::MatrixXd longitudinal(int num_qubits, int site);

Eigen::MatrixXd transverse(int num_qubits, int site);

Eigen::MatrixXd coupling(int num_qubits, int site);

Eigen::MatrixXd hamiltonian(int num_qubits, const Fields& fields);

Eigen::VectorXd magnetization(int num_qubits);
	
#endif /* operators_hpp */

