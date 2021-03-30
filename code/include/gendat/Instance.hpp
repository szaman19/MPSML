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

template <typename T>
class Instance  
{
	public:
		Instance(const Fields<T>& field, const Solver<T>& solver);

		void append_to_file(std::string fpath) const;

		void print(void) const;

		void sign_swap(void);

	private:
		Fields<T> fields;
		Solver<T> wavefx;
};

	
#endif /* Instance_hpp */

