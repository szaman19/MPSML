///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Reader.hpp ***                              //
//                                                                           //
// created December 14, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Reader_hpp
#define Reader_hpp

#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <gendat/Fields.hpp>

template <typename T>
class Reader  
{
	public:
		Reader(int num_qubits, std::string fpath) 
			: num_qubits(num_qubits), fpath(fpath) { }

		void read(int instance);

	private:
		void read(void);
		int num_qubits;
		std::string fpath;
		std::vector<Fields<T>> fields;
		std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> values;
		std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> wavefx;
};



#endif /* Reader_hpp */

