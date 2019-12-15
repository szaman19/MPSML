///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Parser.hpp ***                              //
//                                                                           //
/// Responsible for parsing input file                                       //
//                                                                           //
// created November 15, 2018                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Parser_hpp
#define Parser_hpp

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include "nixio.hpp"

class Parser
{

	public:
		Parser(std::string inputfile);

		bool job(std::string jobname);

		std::string value_of_key(std::string key);

		std::vector<std::string> user_strings(void);

	private:
		std::vector<std::string> file_strings;
		std::multimap<std::string, std::string> inputs; 

		void print(void);
		bool is_comment(std::string key);
		void capture_input(std::istringstream &streamline);
		std::string quote_delimited_field(std::istringstream &streamline);
};

#endif /* Parser_hpp */
