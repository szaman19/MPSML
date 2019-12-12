///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Parser.cpp ***                              //
//                                                                           //
/// Responsible for parsing input file                                       //
//                                                                           //
// created November 15, 2018                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Parser.hpp"

Parser::Parser(std::string inputfile) 
{
	std::fstream file(inputfile.c_str());
	std::string line;
	std::istringstream streamline;

	if (file.is_open()) {
		while (getline(file, line)) file_strings.push_back(line);
		
		for (unsigned i = 0; i < file_strings.size(); ++i) {
			form_streamline(file_strings[i], streamline);
			capture_input(streamline);
		}

	} else {
		open_ascii_escape("red");
		std::cout << "\n ERROR!!! COULD NOT OPEN " << inputfile << "\n" << std::endl;
		close_ascii_escape();
		std::exit(-1);
	}
}

void Parser::capture_input(std::istringstream &streamline)
{
	std::string key;
	std::string eql;
	std::string val;	

	streamline >> key;
	streamline >> eql;

	if (key == "user_string") 
		val = quote_delimited_field(streamline);
 	else 
		streamline >> val;
	
	if (key.size() != 0 && val.size() != 0 && !is_comment(key)) 
		inputs.insert(std::make_pair(key, val));
}

std::string Parser::quote_delimited_field(std::istringstream &streamline)
{
	std::string delimiter = "\"";
	
	std::string tmp;
	std::vector<std::string> us;
	
	streamline >> tmp;

	while (true) {
		streamline >> tmp;
		if (tmp == delimiter) break;
		us.push_back(tmp + " ");
	}

	tmp = "";

	for (auto i : us) tmp += i;

	return tmp;
}

bool Parser::is_comment(std::string key)
{
	std::string comment_specifier = "#";

	if (key.find(comment_specifier) != std::string::npos)
		return true;
	else 
		return false;
}

void Parser::print(void)
{
	auto it = inputs.begin();

	while (it != inputs.end()) {
		std::cout << it->first << " = " << it->second << "\n"; 
		++it;
	}
}

bool Parser::job(std::string jobname)
{
	auto it = inputs.find(jobname);

	if (it != inputs.end() && it->second == "true")
		return true;
	else
		return false;
}

std::vector<std::string> Parser::user_strings(void)
{
	std::vector<std::string> v;
	std::string key = "user_string";

	for (auto it = inputs.begin(); it != inputs.end(); ++it)
	{
		if (it->first == key)
		{
			v.push_back(it->second);
		}
	}

	return v;
}

std::string Parser::value_of_key(std::string key)
{
	auto it = inputs.find(key);

	if (it != inputs.end()) 
		return it->second;
	else
		return "";
}











