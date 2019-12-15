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
#include "Fields.hpp"

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

template <typename T>
void Reader<T>::read(int instance)
{
	int dim = (int)(pow(2, num_qubits) + 0.5);	

	std::streamoff offset = 3*num_qubits + dim + dim*dim;

	std::ifstream file(fpath.c_str(), std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "ERROR: COULD NOT OPEN " << fpath << '\n';
		exit(-1);
	}

	Fields<T> tmp_fields(num_qubits);	
	Eigen::Matrix<T, Eigen::Dynamic, 1> tmp_values(dim);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp_wavefx(dim, dim);
	
	file.seekg(offset * instance * sizeof(T), std::ios::beg);

	if (file.eof())
	{
		std::cerr << "ERROR: FILESTREAM ERROR\n\n";
		exit(-1);
	}

	file.read((char*)tmp_fields.coupling.data(), num_qubits * sizeof(T));
	file.read((char*)tmp_fields.transverse.data(), num_qubits * sizeof(T));
	file.read((char*)tmp_fields.longitudinal.data(), num_qubits * sizeof(T));
	file.read((char*)tmp_values.data(), dim * sizeof(T));
	file.read((char*)tmp_wavefx.data(), dim * dim * sizeof(T));

	fields.push_back(tmp_fields);
	values.push_back(tmp_values);
	wavefx.push_back(tmp_wavefx);

	fields[0].print();
	std::cout << values[0].transpose() << "\n\n";	
	std::cout << wavefx[0] << "\n\n";
}


#endif /* Reader_hpp */

