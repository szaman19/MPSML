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
#include <fstream>
#include <vector>
#include <random>
#include <string>

class Fields
{
	public: 
		Fields() {};

		Fields(int num_qubits, double J, double Bx, double Bz, double W);

		void append_to_file(std::string filename) const;
		void read(std::string fpath, std::streampos pos = std::ios::beg);

		void print(void) const;

		double coupling;
		double transverse;
		double longitudinal;
		double disorder_strength;
		std::vector<double> local_fields;
};

	
#endif /* Fields_hpp */

