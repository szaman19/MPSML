///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Optimizer.hpp ***                             //
//                                                                           //
// created December 09, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Optimizer_hpp
#define Optimizer_hpp

#include <string>
#include <functional>
#include "Network.hpp"

typedef std::function<void(Network& network)> update_fcn_ptr;

class Optimizer  
{
	public:
		Optimizer(std::string optimizer_type,
				double learning_rate = 0.1, 
				double decay_rate = 0.90, 
				double epsilon_conditioner = 1e-6);

		update_fcn_ptr update;

		void set_learning_rate(double value);
		void set_decay_rate(double value);
		void set_epsilon_conditioner(double value);

	private:
		bool m_buffers_initialized = false;
		double m_learning_rate;
		double m_decay_rate; 
		double m_epsilon_conditioner;

		void set_update_pointer(std::string optimizer_type);

		void initialize_accumulation_buffers(const Network& network);
		void compute_gradients(const Network& network);
		void accumulate_gradient(void);
		void compute_update(void);
		void accumulate_update(void);
		void apply_update(Network& network);

		Network m_gradients;
		Network m_updates;
		Network m_gradient_accumulator;
		Network m_update_accumulator;

		void stochastic_gradient_descent(Network& network);
		void adaptive_delta(Network& network);
};

	
#endif /* Optimizer_hpp */

