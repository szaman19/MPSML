///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                         *** Activation.cpp ***                            //
//                                                                           //
// created June 18, 2019                                                     //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Activation.hpp"

Activation::Activation(std::string activation_type)
{
	if (activation_type == "linear")
	{
		activation = linear_0;
		derivative = linear_1;
	}
	else if (activation_type == "tanh")
	{
		activation = tanh_0;
		derivative = tanh_1;
	}
	else if (activation_type == "sigmoid")
	{
		activation = sigmoid_0;
		derivative = sigmoid_1;
	}
	else if (activation_type == "relu")
	{
		activation = relu_0;
		derivative = relu_1;
	}
	else if (activation_type == "elu")
	{
		activation = elu_0;
		derivative = elu_1;
	}
	else
	{
		open_ascii_escape("red");
		std::cerr << "\nERROR!: I DON'T KNOW THAT ACTIVATION FUNCTION\n\n";
		close_ascii_escape();	
		exit(-1);
	}
}

double Activation::relu_0(double x)
{
	if (x < 0) return 0;
	else return x;
}

double Activation::relu_1(double x)
{
	if (x < 0) return 0;
	else return 1;
}

double Activation::elu_0(double x)
{
	if (x < 0) return ELU_ALPHA * (std::exp(x) - 1);
	else return x;
}

double Activation::elu_1(double x)
{
	if (x < 0) return ELU_ALPHA * std::exp(x);
	else return 1;
}

double Activation::sigmoid_0(double x)
{
	return 1.0 / (1.0 + std::exp(x));
}

double Activation::sigmoid_1(double x)
{
	return sigmoid_0(x) * (1.0 - sigmoid_0(x));
}

double Activation::tanh_0(double x)
{
	return std::tanh(x);
	//return x * (10395.0 + 1260.0*x*x + 21.0*x*x*x*x) / 
		//(10395.0 + 4725.0*x*x + 210.0*x*x*x*x + x*x*x*x*x*x);
}

double Activation::tanh_1(double x)
{
	return 1.0 - std::tanh(x) * std::tanh(x);
	//return 1.0 - tanh_0(x) * tanh_0(x);
	//return (-21.0*x*x*x*x*x*x*x*x*x*x  + 
			//630.0*x*x*x*x*x*x*x*x - 
			//18900.0*x*x*x*x*x*x + 
			//496125.0*x*x*x*x - 
			//9823270.0*x*x + 
			//108056000.0) / 
			//( (x*x*x*x*x*x + 210.0*x*x*x*x + 4725.0*x*x + 10395.0) * 
			  //(x*x*x*x*x*x + 210.0*x*x*x*x + 4725.0*x*x + 10395.0) );
}

double Activation::linear_0(double x)
{
	return x;
}

double Activation::linear_1(double)
{
	return 1;
}
