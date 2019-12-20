///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            *** Model.hpp ***                              //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Model_hpp
#define Model_hpp

#include <vector>
#include <phynet/Network.hpp>
#include <phynet/Loss.hpp>
#include <phynet/Optimizer.hpp>
#include <gendat/Operators.hpp>

template <typename T>
class Model  
{
	public:
		Model(const std::vector<Network<T>>& networks, const Loss<T>& loss, 
			  const Optimizer<T>& optimizer);

		void learn_from(const Dataset<T>& dataset); 

		T mse(const Dataset<T>& dataset);

		T predictive_power(const Dataset<T>& dataset, int epoch);

		Eigen::RowVectorXd pure_cost(const Dataset<T>& dataset);

		//void read(std::string fpath);
		//void save(std::string fpath);

		//void write_coefficients(const Dataset<T> &dataset, std::string fpath, int epoch);
		//void write_schrodinger_error(const Dataset<T> &dataset, std::string fpath, int epoch);
		//void write_lyapunov_estimate(const Dataset<T> &dataset, std::string fpath, int epoch);
		//void write_average_magnetization(const Dataset<T> &dataset, std::string fpath, int epoch);

	private:
		std::vector<Network<T>> networks;
		Loss<T> loss; 
		Optimizer<T> optimizer;		

		void pretty_print(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
};
	
#endif /* Model_hpp */

