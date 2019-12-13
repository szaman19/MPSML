///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** pauli.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "pauli.hpp"

Eigen::Matrix2d eye()
{
	Eigen::Matrix2d I;
	I << 1, 0, 0, 1;
	return I;
}

Eigen::Matrix2d sigmax()
{
	Eigen::Matrix2d x;
	x << 0, 1, 1, 0;
	return x;
}

Eigen::Matrix2d sigmaz()
{
	Eigen::Matrix2d z;
	z << 1, 0, 0, -1;
	return z;
}

Eigen::Matrix2cd sigmay()
{
	std::complex<double> i(0,1);
	Eigen::Matrix2cd y;
	y << 0, -i, i, 0;	
	return y;	
}
