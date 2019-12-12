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

#include "Network.hpp"
#include "Loss.hpp"
#include "Optimizer.hpp"

class Model  
{
	public:
		Model(const Network& network, const Loss& loss, const Optimizer& optimizer);

		void learn_from(const Dataset& dataset); 
		double mse(const Dataset& dataset);

		void read(std::string filename);
		void save(std::string filename);

		void write_coefficients(const Dataset &dataset, std::string filename, int epoch);
		void write_schrodinger_error(const Dataset &dataset, std::string filename, int epoch);
		void write_lyapunov_estimate(const Dataset &dataset, std::string filename, int epoch);
		void write_average_magnetization(const Dataset &dataset, std::string filename, int epoch);

	private:
		Network m_net;
		Loss m_loss; 
		Optimizer m_optimizer;		
};

	
#endif /* Model_hpp */

