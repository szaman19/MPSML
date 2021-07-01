///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Reader.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <dump/Reader.hpp>

template <typename T>
void Reader<T>::read(void)
{
	int dim = (int)(pow(2, num_qubits) + 0.5);	

	long offset = (3*num_qubits + 1 + dim)*sizeof(T);
	int num_instances = boost::filesystem::file_size(fpath)/offset;

	std::ifstream file(fpath.c_str(), std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "ERROR: COULD NOT OPEN " << fpath << '\n';
		exit(-1);
	}

	Fields<T> tmp_fields(num_qubits);	
	Eigen::Matrix<T, Eigen::Dynamic, 1> tmp_values(1);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp_wavefx(dim,1);
	
	if (file.eof())
	{
		std::cerr << "ERROR: FILESTREAM ERROR\n\n";
		exit(-1);
	}

	for (int i = 0; i < num_instances; ++i)
	{
		file.read((char*)tmp_fields.coupling.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_fields.transverse.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_fields.longitudinal.data(), num_qubits * sizeof(T));
		file.read((char*)tmp_values.data(), sizeof(T));
		file.read((char*)tmp_wavefx.data(), dim * sizeof(T));

		fields.push_back(tmp_fields);
		values.push_back(tmp_values);
		wavefx.push_back(tmp_wavefx);
	}

}

template <typename T>
void Reader<T>::print(void) const
{
	for (std::size_t i = 0; i < fields.size(); ++i)
	{
		std::cout << std::scientific;
		std::cout << "--------------------------- INSTANCE " << i+1 << " ---------------------------\n";
		fields[i].print();

		std::cout << "Eigenvalues:\n";
		std::cout << values[i].transpose() << "\n\n";
		std::cout << "Eigenvectors:\n";
		std::cout << wavefx[i] << "\n";
		std::cout << "--------------------------------------------------------------------\n\n\n";
	}
}

template class Reader<double>;
