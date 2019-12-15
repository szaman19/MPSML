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
#include "Network.hpp"

typedef std::function<void(Network&, const Dataset&, std::size_t)> loss_fcn_ptr;

class Loss  
{
	public:
		Loss(std::string loss, 
				double lagrange_multiplier = 0.1, 
				double trade_off_parameter = 0.1, 
				double random_domain_bound = 0.1);

		loss_fcn_ptr compute;

		void set_lagrange_multiplier(double value);
		void set_trade_off_parameter(double value);
		void set_random_domain_bound(double value);

	private:
		double m_lagrange_multiplier;
		double m_trade_off_parameter;
		double m_random_domain_bound;

		void set_compute_pointer(std::string loss);

		void quadratic(Network& n, const Dataset& d, std::size_t batch);
		void quadratic_plus_schrodinger(Network& n, const Dataset& d, std::size_t batch);
		void physics_perturbed_quadratic(Network& n, const Dataset& d, std::size_t batch);
		void randomly_perturbed_quadratic(Network& n, const Dataset& d, std::size_t batch);
};

	
#endif /* Loss_hpp */

