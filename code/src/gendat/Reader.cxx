///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Reader.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <gendat/Reader.hpp>

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
}

template <typename T>
void Reader<T>::print(void) const
{

	for (std::size_t i = 0; i < fields.size(); ++i)
	{
		fields[i].print();

		std::cout << "Eigenvalues:\n";
		std::cout << values[i].transpose() << "\n\n";
		std::cout << "Eigenvectors:\n";
		std::cout << wavefx[i] << "\n\n";
	}
}

template class Reader<double>;
