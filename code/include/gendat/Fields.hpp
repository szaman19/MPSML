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


#endif /* Fields_hpp */

