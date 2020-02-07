///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                           *** Loss.cpp ***                                //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Loss_hpp
#define Loss_hpp

#include <cmath>
#include <string>
#include <functional>
#include <phynet/Network.hpp>
#include <gendat/Operators.hpp>

template <typename T>
using Batch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using NetVec = std::vector<Network<T>>;

template <typename T>
class Loss  
{
	typedef std::function<void(NetVec<T>&, const Dataset<T>&, int batch)> loss_fcn_ptr;
	
	public:
		Loss(std::string loss, const Operators<T>& operators,
				T lagrange_multiplier = 0.1, 
				T trade_off_parameter = 0.1, 
				T random_domain_bound = 0.1);

		loss_fcn_ptr compute;

		Eigen::RowVectorXd pure_cost(NetVec<T>& nets, const Dataset<T>& data);

		void set_lagrange_multiplier(T value);
		void set_trade_off_parameter(T value);
		void set_random_domain_bound(T value);

	private:
		T lagrange_multiplier;
		T trade_off_parameter;
		T random_domain_bound;

		Operators<T> operators;
		
		Eigen::SparseMatrix<T> lagrange_matrix;

		void set_compute_pointer(std::string loss);

		void quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void quadratic_plus_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void unsupervised_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void supervised_schrodinger(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void abs_formulation(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void physics_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void randomly_perturbed_quadratic(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void rayleigh_ritz(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void shuijun_diannong(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void data_shuijun_diannong(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void raw_sums(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void sigmoid_unitarity(NetVec<T>& nets, const Dataset<T>& data, int batch);
		void sigmoid_frobenius(NetVec<T>& nets, const Dataset<T>& data, int batch);
};

inline double sigmoid(double x)
{
	return 1.0 / (1.0 + std::exp(-x));
}
	
#endif /* Loss_hpp */

