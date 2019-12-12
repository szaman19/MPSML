#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

const int NUM_FIELDS = 12000;
const int NUM_EPOCHS = 1000;

double bb[NUM_FIELDS][NUM_EPOCHS];
double pg[NUM_FIELDS][NUM_EPOCHS];
double c2[NUM_FIELDS][NUM_EPOCHS];
double m0[NUM_FIELDS];
double bx[NUM_FIELDS];
double avg_over_field[NUM_EPOCHS];
double avg_over_epoch[NUM_FIELDS];

string base_dir = "../../data/output/J-equal-minus-one/single-phase/";
string num_site = "4-site-";
string in_out = "-hamil-gs-coeff-epoch=";
string tail = ".dat";


inline void form_streamline(std::string &line, std::istringstream &streamline);
void get_predicted_mz(string algo, double data[NUM_FIELDS][NUM_EPOCHS]);
void get_actual_mz(void);
void write_avg_over_epoch(string fpath, double data[NUM_FIELDS][NUM_EPOCHS]);
void write_avg_over_field(string fpath, double data[NUM_FIELDS][NUM_EPOCHS]);

int main( /* int argc, char *argv[] */ )
{
	get_actual_mz();

	get_predicted_mz("bb", bb);
	get_predicted_mz("pg", pg);
	get_predicted_mz("c2", c2);

	write_avg_over_epoch(base_dir + "bb-avg-over-epoch.txt", bb);
	write_avg_over_epoch(base_dir + "pg-avg-over-epoch.txt", pg);
	write_avg_over_epoch(base_dir + "c2-avg-over-epoch.txt", c2);

	write_avg_over_field(base_dir + "bb-avg-over-field.txt", bb);
	write_avg_over_field(base_dir + "pg-avg-over-field.txt", pg);
	write_avg_over_field(base_dir + "c2-avg-over-field.txt", c2);

	return 0;
}

void get_actual_mz(void)
{
	ifstream file(base_dir + "4-site-Mz.txt");
	string line;
	istringstream streamline;

	int c = 0;

	while (getline(file, line))
	{
		form_streamline(line, streamline);
		streamline >> bx[c];
		streamline >> m0[c];
		++c;
	}

}

void write_avg_over_field(string fpath, double data[NUM_FIELDS][NUM_EPOCHS])
{
	double tmp;
	ofstream file(fpath);

	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
	{
		tmp = 0;
		for (int field = 0; field < NUM_FIELDS; ++field)
		{
			tmp += sqrt((m0[field] - data[field][epoch])*(m0[field] - data[field][epoch]));
		}

		file << scientific << setw(20);
		file << epoch << '\t';
		file << 1 - tmp / (2 * NUM_FIELDS) << '\n';
	}
}

void write_avg_over_epoch(string fpath, double data[NUM_FIELDS][NUM_EPOCHS])
{
	double tmp;

	ofstream file(fpath);

	for (int field = 0; field < NUM_FIELDS; ++field)
	{
		tmp = 0;
		for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
		{
			tmp += sqrt((m0[field] - data[field][epoch])*(m0[field] - data[field][epoch]));
		}

		file << scientific << setw(20);
		file << bx[field] << '\t';
		file << 1 - tmp / (2 * NUM_EPOCHS) << '\n';
	}
}


void get_predicted_mz(string algo, double data[NUM_FIELDS][NUM_EPOCHS])
{
	string fpath;
	ifstream file;
	istringstream streamline;
	string line;
	double bx, bz, c1, cn;
	int count;

	for (int epoch = 0; epoch < NUM_EPOCHS; ++epoch)
	{
		fpath = base_dir + num_site + algo + in_out + to_string(epoch) + tail;
		file.open(fpath);
		if (file.is_open())
		{
			getline(file, line);	
			count = 0;
			while (getline(file, line))
			{
				form_streamline(line, streamline);	
				streamline >> bx;
				streamline >> bz;
				streamline >> c1;

				for (int i = 0; i < 15; ++i) streamline >> cn;
			
				data[count++][epoch] = c1*c1 - cn*cn;
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

