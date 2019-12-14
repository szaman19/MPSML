///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Instance.hpp ***                              //
//                                                                           //
// created December 13, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Instance_hpp
#define Instance_hpp

#include <Eigen/Dense>
#include "Fields.hpp"
#include "Solver.hpp"
#include <boost/filesystem.hpp>

class Instance  
{
	public:
		Instance() {};
		Instance(const Fields& field, const Solver& solver);

		const Eigen::VectorXd& spectrum(void) const;
		const Eigen::MatrixXd& unitary(void) const;

		double eigenenergy(int state_index) const;

		void append_to_file(std::string fpath) const;

		void print(void) const;

		void read(std::string fpath, std::streampos pos = std::ios::beg);

	private:
		Fields fields;
		Solver wavefx;
};


	
#endif /* Instance_hpp */

