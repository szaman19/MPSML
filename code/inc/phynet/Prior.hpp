///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Prior.hpp ***                               //
//                                                                           //
// created June 30, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Prior_hpp
#define Prior_hpp

#include "Eigen/Dense"
#include "nixio.hpp"
#include <vector>

typedef Eigen::Map<const Eigen::ArrayXd> MaptoArray;
typedef Eigen::Map<const Eigen::MatrixXd> MaptoFMat;

class Prior  
{
	public:
		Prior(std::string inputfile, double intercessor_value); 

		const Eigen::MatrixXd& operant(std::size_t instance) const;
		const Eigen::MatrixXd& eigenvalue(std::size_t instance) const;
		const Eigen::MatrixXd& intercessor(std::size_t instance) const;

		std::size_t inner(void)     const { return m_inner_dimension; }
		std::size_t outer(void)     const { return m_outer_dimension; }
		std::size_t instances(void) const { return m_instances; }

	private:
		std::size_t m_instances;
		std::size_t m_inner_dimension;
		std::size_t m_outer_dimension;
		
		std::vector<Eigen::MatrixXd> m_operants;
		std::vector<Eigen::MatrixXd> m_eigenvalues;
		std::vector<Eigen::MatrixXd> m_intercessor;
};

	
#endif /* Prior_hpp */

