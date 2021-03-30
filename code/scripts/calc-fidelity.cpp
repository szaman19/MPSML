#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

const int NUM_FIELDS = 12000;
const int NUM_EPOCHS = 1001;
const int NUM_COEFFS = 16 + 2; // include two field values 

double data[NUM_EPOCHS][NUM_FIELDS][NUM_COEFFS];

string base_dir = "../../data/output/J-equal-minus-one/single-phase/";
string num_site = "4-site-";
string in_out = "-hamil-gs-coeff-epoch=";
string tail = ".dat";

struct Row { vector<double> coeffs; };
bool compare_row(const Row &row1, const Row &row2) ;
inline void form_streamline(std::string &line, std::istringstream &streamline);
void get_actual_mz(void);
void get_data(string algo);
void sort_single_epoch(int epoch);
void sort_all_epochs(void);
double fidelity(int epoch, int field);
double susceptibility(int epoch, int field);
void write_susceptibility(string algo);
void norm_wavefunctions(void);

vector<string> algos = {"bb", "pg", "c2"};

int main( /* int argc, char *argv[] */ )
{

	for (auto algo : algos)
	{
		get_data(algo);
		sort_all_epochs();
		norm_wavefunctions();
		write_susceptibility(algo);
	}

	return 0;
}

void norm_wavefunctions(void)
{
	double tmp;
	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
	{
		for (int field = 0; field < NUM_FIELDS; ++field)
		{
			tmp = 0;
			for (int coeff = 2; coeff < NUM_COEFFS; ++coeff)
				tmp += data[epoch][field][coeff] * data[epoch][field][coeff];

			for (int coeff = 2; coeff < NUM_COEFFS; ++coeff)
				data[epoch][field][coeff] /= sqrt(tmp);
		}
	}
}

void write_susceptibility(string algo)
{
	ofstream file;
	string fpath;

	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
	{
		fpath = base_dir + num_site + algo + "-hamil-gs-fidelity-epoch=" + to_string(epoch) + tail;
		file.open(fpath);
		
		file << scientific;
		file << setw(10);

		for (int field = 0; field < NUM_FIELDS - 1; ++field)
		{
			file << data[epoch][field][0] << '\t';
			file << susceptibility(epoch, field) << '\n';
		}

		file.close();
	}
}

double susceptibility(int epoch, int field)
{
	if (field == NUM_FIELDS || epoch == NUM_EPOCHS)
	{
		cerr << "out of range...\n" << endl;
		exit(-10);
	}
	double denom = data[epoch][field][0] - data[epoch][field + 1][0];

	denom = denom * denom;

	return (2 - 2 * fidelity(epoch, field)) / denom;
}


double fidelity(int epoch, int field)
{
	double tmp = 0;
		
	for (int coeff = 2; coeff < NUM_COEFFS; ++coeff)
		tmp += data[epoch][field][coeff] * data[epoch][field + 1][coeff];

	return abs(tmp);
}


void sort_all_epochs(void)
{
	cout << "sorting...\n";
	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
		sort_single_epoch(epoch);
}

void sort_single_epoch(int epoch)
{
	Row row;
	vector<Row> rows;

	for (int field = 0; field < NUM_FIELDS; ++field)
	{
		for (int coeff = 0; coeff < NUM_COEFFS; ++coeff)
		{
			row.coeffs.push_back(data[epoch][field][coeff]);
		}
		rows.push_back(row);
		row.coeffs.clear();
	}

	sort(&rows[0], &rows[NUM_FIELDS], compare_row);

	for (int field = 0; field < NUM_FIELDS; ++field)
	{
		for (int coeff = 0; coeff < NUM_COEFFS; ++coeff)
		{
			data[epoch][field][coeff] = rows[field].coeffs[coeff];
		}
	}
}

void get_data(string algo)
{
	cout << "reading...\n";
	string fpath, line;
	ifstream file;
	istringstream streamline;

	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
	{
		fpath = base_dir + num_site + algo + in_out + to_string(epoch) + tail;
		file.open(fpath);
		if (file.is_open())
		{
			getline(file, line);	

			for (int field = 0; field < NUM_FIELDS; ++field)
			{
				getline(file, line);
				form_streamline(line, streamline);

				for (int coeff = 0; coeff < NUM_COEFFS; ++coeff)
					streamline >> data[epoch][field][coeff];
			}	
		}
		file.close();
	}
}

inline void form_streamline(std::string &line, std::istringstream &streamline)
{
	streamline.str(std::string());
	streamline.clear();
	streamline.str(line);
}

bool compare_row(const Row &row1, const Row &row2) 
{
	return row1.coeffs[0] < row2.coeffs[0];
}
