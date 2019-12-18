///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                      *** Instance.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <gendat/Instance.hpp>

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

template class Instance<double>;
