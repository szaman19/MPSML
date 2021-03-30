#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

int main( /* int argc, char *argv[] */ )
{
	cout << sizeof(double) << '\n';
	cout << sizeof(float) << '\n';

	cout << boost::filesystem::file_size("4-qubits.bin") << '\n';


	return 0;
}
