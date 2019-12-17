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
#include <phynet/Network.hpp>

template <typename T>
class Optimizer  
{
	typedef std::function<void(Network<T>& network)> update_fcn_ptr;

	public:
		Optimizer(std::string optimizer_type,
				T learning_rate = 0.1, 
				T decay_rate = 0.90, 
				T epsilon_conditioner = 1e-6);

		update_fcn_ptr update;

		void set_learning_rate(T value);
		void set_decay_rate(T value);
		void set_epsilon_conditioner(T value);

	private:
		bool m_buffers_initialized = false;
		T m_learning_rate;
		T m_decay_rate; 
		T m_epsilon_conditioner;

		void set_update_pointer(std::string optimizer_type);

		void initialize_accumulation_buffers(const Network<T>& network);
		void compute_gradients(const Network<T>& network);
		void accumulate_gradient(void);
		void compute_update(void);
		void accumulate_update(void);
		void apply_update(Network<T>& network);

		Network<T> m_gradients;
		Network<T> m_updates;
		Network<T> m_gradient_accumulator;
		Network<T> m_update_accumulator;

		void stochastic_gradient_descent(Network<T>& network);
		void adaptive_delta(Network<T>& network);
};

	
#endif /* Optimizer_hpp */

