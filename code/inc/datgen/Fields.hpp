///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Fields.hpp ***                              //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Fields_hpp
#define Fields_hpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <string>

template <typename T>
class Fields
{
	public: 
		Fields() {};
		Fields(int num_qubits);
		Fields(int num_qubits, T J, T Bx, T Bz, T W);

		void append_to_file(std::string filename) const;
		void read(std::string fpath, std::streampos pos = std::ios::beg);

		void print(void) const;

		std::vector<T> coupling;
		std::vector<T> transverse;
		std::vector<T> longitudinal;
};

template <typename T>
Fields<T>::Fields(int num_qubits, T J, T Bx, T Bz, T W)
{
	std::random_device rd;
	std::default_random_engine engine(rd());
	std::uniform_real_distribution<T> dist(-W, W);

	for (int i = 0; i < num_qubits; ++i) 
	{
		coupling.push_back(J + dist(engine));
		transverse.push_back(Bx + dist(engine));
		longitudinal.push_back(Bz + dist(engine));
	}
}

template <typename T>
void Fields<T>::print(void) const
{
	std::cout << '\n';
	for (int i = 0; i < coupling.size(); ++i)
	{
		std::cout << " J[" << i << "] = " << std::setw(10) << coupling[i] << '\t';
		std::cout << "Bx[" << i << "] = " << std::setw(10) << transverse[i] << '\t';
		std::cout << "Bz[" << i << "] = " << std::setw(10) << longitudinal[i] << '\n';
	}
	std::cout << '\n';
}

template <typename T>
Fields<T>::Fields(int num_qubits)
{
	coupling.resize(num_qubits);
	transverse.resize(num_qubits);
	longitudinal.resize(num_qubits);
}

template <typename T>
void Fields<T>::append_to_file(std::string filename) const
{
	std::ofstream file(filename.c_str(), std::ios::app | std::ios::binary);	

	file.write((char*)coupling.data(), coupling.size() * sizeof(T));
	file.write((char*)transverse.data(), transverse.size() * sizeof(T));
	file.write((char*)longitudinal.data(), longitudinal.size() * sizeof(T));
}
	
#endif /* Fields_hpp */

