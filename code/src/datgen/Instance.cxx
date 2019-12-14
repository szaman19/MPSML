///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                      *** Instance.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Instance.hpp"

Instance::Instance(const Fields& field, const Solver& solver) 
	: fields(field), wavefx(solver)
{
	// null body
}

const Eigen::MatrixXd& Instance::unitary(void) const
{
	return wavefx.eigenvectors();
}

const Eigen::VectorXd& Instance::spectrum(void) const
{
	return wavefx.eigenvalues();
}

double Instance::eigenenergy(int state_index) const 
{
	return wavefx.eigenvalues()(state_index);
}

void Instance::append_to_file(std::string fpath) const
{
	fields.append_to_file(fpath);	
	wavefx.append_to_file(fpath);
}

void Instance::read(std::string fpath, std::streampos pos)
{
	std::streamoff off = (4 + fields.local_fields.size()) * sizeof(double); 

	fields.read(fpath, pos);
	wavefx.read(fpath, pos + off);
}

void Instance::print(void) const
{
	fields.print();
	wavefx.print();
}
