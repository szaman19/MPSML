///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** datgen main.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <filesystem>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Fields.hpp>
#include <Generator.hpp>
#include <Operators.hpp>

struct GeneratorConfig
{

	GeneratorConfig(){}; 

	int dBx, dBz, qubits;
	double Bx_min, Bz_min, Bx_max, Bz_max, coupling;
	std::string fpath;
 
  friend std::ostream& operator<<(std::ostream &os, const GeneratorConfig& G);

};


constexpr unsigned int str2int(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}



std::ostream& operator<<(std::ostream &os, const GeneratorConfig& G){
	os << "Generator configuration: \n"
		 << "Number of qubits \t" << G.qubits << '\n'
		 << "Number of Bx Partitions \t" << G.dBx << '\n'
		 << "Number of Bz Partitions \t" << G.dBz << '\n'
		 << "Bx minimum \t" << G.Bx_min << '\n'
		 << "Bx maximum \t" << G.Bx_max << '\n'
		 << "Bz minimum \t" << G.Bz_min << '\n'
		 << "Bz maximum \t" << G.Bz_max << '\n'
		 << "Coupling strengh \t" << G.coupling << '\n'
		 << "Output locatoin \t" << G.fpath;
		 return os;
}

void config_parser(GeneratorConfig &configObj){

	std::ifstream config_file; 
	config_file.open("config.ini", std::ios::in);

	assert(("Couldn't open file", config_file.is_open()));
	std::string line; 
	while(getline(config_file, line)){
		

		if (line != "[DEFAULT]"){ // This is the  default genereated by python configparser

			std::string delimiter = "=";

			auto delim_pos = line.find(delimiter);

			std::string arg = line.substr(0, delim_pos-1);
			std::string val = line.substr(delim_pos + 2 , line.length());

			switch(str2int(arg.c_str())){
				case str2int("Qubits"):
					configObj.qubits = std::stoi(val);
					break; 
				case str2int("BxPartitions"):
						configObj.dBx= std::stoi(val);
					break;
				case str2int("BzPartitions"):
						configObj.dBz= std::stoi(val);
					break;
				case str2int("BxMin"):
						configObj.Bx_min= std::stod(val);
					break;
				case str2int("BzMin"):
						configObj.Bz_min= std::stod(val);
					break;
				case str2int("BxMax"):
						configObj.Bx_max = std::stod(val);
					break;
				case str2int("BzMax"):
						configObj.Bz_max = std::stod(val);
					break;
				case str2int("Coupling"):
						configObj.coupling = std::stod(val);
					break;
				case str2int("OutputLoc"):
						// Quick little hack but requires C++17
						if (val == "."){
							val = std::filesystem::current_path().string() + "/out.bin";
						}
						configObj.fpath= val;
					break;
				case str2int("FullSpectrum"):
						std::cout << "For future Shehtab to worry about";
					break;
				default:	
					throw std::invalid_argument(arg + " Invalid argument in config file");
			}

		}
	}
	config_file.close();
	

}



int main(int argc, char *argv[])
{
	GeneratorConfig gen_configs;
	config_parser(gen_configs);

	std::cout << gen_configs << std::endl; 

	
	std::string model = "ising";

	Generator<double>generator;	

	// // required inputs 
	generator.qubits = gen_configs.qubits;
	generator.fpath =  gen_configs.fpath;

	generator.coupling = gen_configs.coupling;
	
	generator.dBx = gen_configs.dBx;
	generator.dBz =  gen_configs.dBx;

	generator.Bx_min =  gen_configs.Bx_min;
	generator.Bx_max = gen_configs.Bx_max;
	
	generator.Bz_min = gen_configs.Bz_min;
	generator.Bz_max = gen_configs.Bz_max;

	generator.run();

	return 0;
}


