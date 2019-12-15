///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** nixio.cpp ***                               //
//                                                                           //
/// Need set of basic nix sys commands bc std::filesystem not portably yet   //
//                                                                           //
// created November 24, 2018                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <phynet/nixio.hpp>

void get_streamline(std::ifstream &file, std::istringstream &streamline)
{
	std::string line;

	std::getline(file, line);
	streamline.str(std::string());
	streamline.clear();
	streamline.str(line);
}

void form_streamline(std::string &line, std::istringstream &streamline)
{
	streamline.str(std::string());
	streamline.clear();
	streamline.str(line);
}

void open_ascii_escape(std::string color)
{
	if (color == "red")    std::cout << "\033[1;31m";
	if (color == "green")  std::cout << "\033[1;32m";
	if (color == "yellow") std::cout << "\033[1;33m";
	if (color == "blue")   std::cout << "\033[1;34m";
	if (color == "purple") std::cout << "\033[1;35m";
}

void close_ascii_escape(void)
{
	std::cout << "\033[0m";
}

void welcome_message()
{
	std::cout << '\n';
	std::cout << "\t------------------------------------------";
	std::cout << '\n';

	std::cout << "\t Starting Physics Guided Machine Learning\n";
	std::cout << "\t          PR_ XX, XXXXXX (XXXX)\n";
	std::cout << "\t------------------------------------------";
	std::cout << '\n';
	std::cout << '\n';
}

int num_lines_in_file(std::string inputfile)
{
	int c = 0;
	std::fstream file(inputfile.c_str());
	std::string line;

	while (getline(file, line)) c++;

	return c;
}
