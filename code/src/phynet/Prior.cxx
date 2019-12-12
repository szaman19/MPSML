///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Prior.cpp ***                               //
//                                                                           //
// created June 30, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Prior.hpp"

Prior::Prior(std::string inputfile, double intercessor_value)
{
	std::ifstream file(inputfile.c_str(), std::ios::binary);

	char *operants_ptr, *eigenvalues_ptr;
	std::size_t operants_bytes, eigenvalues_bytes;

	std::vector<int> header(3);

	std::cout << "Reading prior from " << inputfile << '\n';
	
	if (file.is_open())
	{
		file.read(reinterpret_cast<char*>(header.data()), static_cast<long>(3*sizeof(int)));
		
		m_instances       = static_cast<std::size_t>(header[0]);
		m_inner_dimension = static_cast<std::size_t>(header[1]);
		m_outer_dimension = static_cast<std::size_t>(header[2]);

		if (!"verbose") 
		{
			std::cout << "Got " << m_instances << " instances \n";
			std::cout << "Got " << m_inner_dimension << " for inner dim \n";
			std::cout << "Got " << m_outer_dimension << " for outer dim \n";
		}
		
		//if (m_inner_dimension != m_outer_dimension)
		//{
			//open_ascii_escape("yellow");
			//std::cerr << "\nWARNING!: PRIOR INNER DIMENSION != PRIOR OUTER DIMENSION\n";
			//std::cerr << "I HOPE YOU KNOW WHAT YOU ARE DOING!" << std::endl;
			//close_ascii_escape();
			////exit(-1);
		//}

		Eigen::MatrixXd operant_tmp(m_inner_dimension, m_inner_dimension);
		Eigen::MatrixXd eigenvalue_tmp(m_inner_dimension, m_outer_dimension);
		Eigen::MatrixXd intercessor_tmp(m_inner_dimension, m_outer_dimension);

		// set the diagonals of the intercessor to a constant
		for (long i = 0; i < static_cast<long>(m_inner_dimension); ++i)
		{
			for (long j = 0; j < static_cast<long>(m_outer_dimension); ++j)
			{
				if (i == j) 
				{
					intercessor_tmp(i,i) = intercessor_value;
				}
			}
		}

		// reserve memory for each matrix at each instance 
		for (std::size_t i = 0; i < m_instances; ++i)
		{
			m_operants.push_back(operant_tmp);
			m_eigenvalues.push_back(eigenvalue_tmp);
			m_intercessor.push_back(intercessor_tmp);
		}

		// calculate the number of bytes in each matrix depending on size of prior
		operants_bytes    = sizeof(double) * m_inner_dimension * m_inner_dimension;
		eigenvalues_bytes = sizeof(double) * m_inner_dimension * m_outer_dimension;

		// pull hamiltonians from file 
		for (std::size_t i = 0; i < m_instances; ++i)
		{
			operants_ptr = reinterpret_cast<char*>(m_operants[i].data());
			file.read(operants_ptr, static_cast<long>(operants_bytes));
			if (!file)
			{
				std::cerr << "\nERROR! BROKEN FILE STREAM IN PRIOR\n";
				std::cerr << std::endl;
				exit(-1);
			}
		}

		// pull eigenvalue matrices from file
		for (std::size_t i = 0; i < m_instances; ++i)
		{
			eigenvalues_ptr = reinterpret_cast<char*>(m_eigenvalues[i].data());
			file.read(eigenvalues_ptr, static_cast<long>(eigenvalues_bytes));
			if (!file)
			{
				std::cerr << "\nERROR! BROKEN FILE STREAM IN PRIOR\n";
				std::cerr << std::endl;
				exit(-1);
			}
		}
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: CANNOT OPEN "<< inputfile << '\n' << std::endl;
		close_ascii_escape();
		exit(-1);
	}
}

const Eigen::MatrixXd& Prior::operant(std::size_t instance) const
{
	return m_operants[instance];
}
const Eigen::MatrixXd& Prior::eigenvalue(std::size_t instance) const
{
	return m_eigenvalues[instance];
}
const Eigen::MatrixXd& Prior::intercessor(std::size_t instance) const
{
	return m_intercessor[instance];
}
