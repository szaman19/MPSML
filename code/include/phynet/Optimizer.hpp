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
	typedef std::function<double(void)> learning_schedule_ptr;

	public:
		Optimizer(std::string optimizer_type, std::string learning_schedule_type);

		update_fcn_ptr update;
		learning_schedule_ptr learning_rate;

		void set_learning_rate_min(double value);
		void set_learning_rate_max(double value);
		void set_decay_rate(double value);
		void set_epsilon_conditioner(double value);

		void set_epoch(int epoch);
		void set_epochs(int epochs);
		void set_epochs_per_cycle(int cycles_per_epoch);

	private:
		bool m_buffers_initialized = false;
		T m_learning_rate_min;
		T m_learning_rate_max;
		T m_decay_rate; 
		T m_epsilon_conditioner;

		int m_epoch, m_epochs_per_cycle, m_epochs;

		double constant_schedule();
		double sawtooth_schedule();
		double cosine_schedule();
		double sawtooth_annealed_schedule();
		double cosine_annealed_schedule();

		void set_update_pointer(std::string optimizer_type);
		void set_learning_schedule_pointer(std::string schedule_type);

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

