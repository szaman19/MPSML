///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Fields.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Fields.hpp"

Fields::Fields(int num_qubits, double J, double Bx, double Bz, double W)
	: coupling(J), transverse(Bx), longitudinal(Bz), disorder_strength(W) 
{
	std::random_device rd;
	std::default_random_engine engine(rd());
	std::uniform_real_distribution<double> dist(-W, W);

	for (int i = 0; i < num_qubits; ++i) local_fields.push_back(dist(engine));
	//for (int i = 0; i < num_qubits; ++i) std::cout << local_fields[i] << '\t';
	//std::cout << '\n';
}

void Fields::append_to_file(std::string filename) const
{
	std::ofstream file(filename.c_str(), std::ios::app | std::ios::binary);	

	file.write((char*)&coupling, sizeof(double));
	file.write((char*)&transverse, sizeof(double));
	file.write((char*)&longitudinal, sizeof(double));
	file.write((char*)&disorder_strength, sizeof(double));
	file.write((char*)local_fields.data(), local_fields.size() * sizeof(double));

	//file << coupling << '\t' 
		 //<< transverse << '\t'
		 //<< longitudinal << '\t'
		 //<< disorder_strength << '\t';

	//for (const auto& i : local_fields) file << i << '\t';

	//file << '\n';
}

void Fields::read(std::string fpath, std::streampos pos)
{
	std::ifstream file(fpath.c_str(), std::ios::binary);

	if (file.is_open())
	{
		file.seekg(pos);
		file.read((char*)&coupling, sizeof(double));
		file.read((char*)&transverse, sizeof(double));
		file.read((char*)&longitudinal, sizeof(double));
		file.read((char*)&disorder_strength, sizeof(double));
		file.read((char*)local_fields.data(), local_fields.size()*sizeof(double));
	}
	else
	{
		std::cerr << "ERROR: FAILED TO OPEN " << fpath << '\n';
		exit(-1);
	}
}

void Fields::print(void) const
{
	std::cout << "J:  " << coupling << '\t'
			  << "Bx: " << transverse << '\t'
			  << "Bz: " << longitudinal << '\t'
			  << "W:  " << disorder_strength << '\t';

	std::cout << "hi: " << '\t';
	for (auto i : local_fields) std::cout << i << '\t';
	std::cout << "\n\n";
}
