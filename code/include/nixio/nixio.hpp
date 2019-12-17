///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** nixio.hpp ***                               //
//                                                                           //
/// Conversion of strings and string streams into more useful types          //
//                                                                           //
// created November 24, 2018                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef nixio_hpp
#define nixio_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

/// Pull a stringstream directly from a file by line
void get_streamline(std::ifstream &file, std::istringstream &streamline);

/// Form a stringstream given the corresponding string
void form_streamline(std::string &line, std::istringstream &streamline);

void welcome_message();

void open_ascii_escape(std::string color);
void close_ascii_escape();


int num_lines_in_file(std::string inputfile);
	
#endif /* nixio_hpp */
