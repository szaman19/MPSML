///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                        *** Decomp.cxx ***                                 //
//                                                                           //
// created December 12, 2019                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Decomp.hpp"

Decomp::Decomp(int num_qubits, int num_transverse_fields, int num_disorder_realizations)
{
	Eigen::ArrayXf bx_values = 
		Eigen::ArrayXf::LinSpaced(num_transverse_fields, 0.1, 2);

	double W;

	if (num_disorder_realizations > 1) 
		W = 0.01;
	else
		W = 0.0;

	saes_t waves;

	Fields field;

	for (int iBx = 0; iBx < bx_values.size(); ++iBx)
	{
		for (int iW = 0; iW < num_disorder_realizations; ++iW)
		{
			field = Fields(num_qubits, bx_values[iBx], W);
			waves.compute(hamiltonian(num_qubits, field));

			if (waves.info() != Eigen::Success)
			{
				std::cerr << "failed to diagonalize\n";
				exit(-1);
			}

			wavefx.push_back(waves);
			fields.push_back(field);
		}
	}
}





