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

	private:
		Fields<T> fields;
		Solver<T> wavefx;
};

template <typename T>
Instance<T>::Instance(const Fields<T>& field, const Solver<T>& solver) 
	: fields(field), wavefx(solver)
{
	// null body
}

template <typename T>
void Instance<T>::append_to_file(std::string fpath) const
{
	fields.append_to_file(fpath);	
	wavefx.append_to_file(fpath);
}

template <typename T>
void Instance<T>::print(void) const
{
	fields.print();
	wavefx.print();
}

	
#endif /* Instance_hpp */

