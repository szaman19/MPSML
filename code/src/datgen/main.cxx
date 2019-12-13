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
#include <string>
#include "Decomp.hpp"

int main( /* int argc, char *argv[] */ )
{
	Decomp d(4, 100, 10);	

	std::cout << sizeof(d) << '\n';	

	return 0;
}


