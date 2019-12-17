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

template <typename T>
using Batch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
class Loss  
{
	typedef std::function<void(Network<T>&, const Batch<T>& target_batch)> loss_fcn_ptr;
	
	public:
		Loss(std::string loss, 
				T lagrange_multiplier = 0.1, 
				T trade_off_parameter = 0.1, 
				T random_domain_bound = 0.1);

		loss_fcn_ptr compute;

		void set_lagrange_multiplier(T value);
		void set_trade_off_parameter(T value);
		void set_random_domain_bound(T value);

	private:
		T m_lagrange_multiplier;
		T m_trade_off_parameter;
		T m_random_domain_bound;

		void set_compute_pointer(std::string loss);

		void quadratic(Network<T>& n, const Batch<T>& target_batch);
		//void quadratic_plus_schrodinger(Network<T>& n, const Batch<T>& target_batch);
		//void physics_perturbed_quadratic(Network<T>& n, const Batch<T>& target_batch);
		//void randomly_perturbed_quadratic(Network<T>& n, const Batch<T>& target_batch);
};
	
#endif /* Loss_hpp */

