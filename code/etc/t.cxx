///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** datgen main.cxx ***                           //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;

//typedef std::size_t count_t;
typedef int count_t;

MatrixXd foo(count_t num)
{
	vector<Matrix2d> eyes(num);

	for (count_t i = 0; i < num; ++i) 
		eyes[i] = Matrix2d::Identity();

	MatrixXd tmp = eyes.back();

	for (count_t i = eyes.size() - 2; i >= 0; --i)
		tmp = kroneckerProduct(eyes[i], tmp).eval();

	return tmp;
}

int main( /* int argc, char *argv[] */ )
{

	cout << foo(2) << '\n';


	return 0;
}

