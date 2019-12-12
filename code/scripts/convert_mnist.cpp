#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <vector>

using namespace std;

void form_streamline(std::string &line, std::istringstream &streamline)
{
	streamline.str(std::string());
	streamline.clear();
	streamline.str(line);
}

int main( /* int argc, char *argv[] */ )
{
	std::string inputfile = "data/input/mnist_test.csv";
	std::string outputfile = "data/input/mnist-test.bin";

	std::ifstream file(inputfile.c_str());
	std::string line;
	std::istringstream streamline;

	double pixel;
	int label;
	char comma;

	Eigen::VectorXd h(10);
	h.setZero();

	std::vector<double> features, targets;

	while (getline(file, line))
	{
		form_streamline(line, streamline);	

		streamline >> label;
		streamline >> comma;

		for (int i = 0; i < 10; ++i) if (label == i) h(i) = 1;
		for (int i = 0; i < 10; ++i) targets.push_back(h(i));

		h.setZero();

		while (streamline)
		{
			streamline >> pixel;
			streamline >> comma;
			features.push_back(pixel);
		}
	}
	std::ofstream of(outputfile.c_str(), std::ios::binary);

	int ii = 10000;
	int ff = 28*28;
	int tt = 10;

	of.write(reinterpret_cast<char*>(&ii), static_cast<long>(sizeof(int)));
	of.write(reinterpret_cast<char*>(&ff), static_cast<long>(sizeof(int)));
	of.write(reinterpret_cast<char*>(&tt), static_cast<long>(sizeof(int)));
	of.write(reinterpret_cast<char*>(features.data()), static_cast<long>(sizeof(double)*features.size()));
	of.write(reinterpret_cast<char*>(targets.data()), static_cast<long>(sizeof(double)*targets.size()));

	return 0;
}
