#include "Utilities.hpp"

#include <boost/math/special_functions.hpp>

namespace Utilities
{

	const double pi = acos(-1.0);

}

double Utilities::gaussianRand()
{
	double result;
	do
	{
		double u1 = (float)rand()/(float)RAND_MAX;
		double u2 = (float)rand()/(float)RAND_MAX;
		result = std::sqrt((double)-2.0*std::log(u1))*std::cos(2.0*pi*u2);
	}
	while ( boost::math::isnan((double)result)==true || boost::math::isinf(result)==true  || result >1 || result <-1);
	
	if (result != result)
	{
		std::cout << "ERR: gaussian rand result: " << "(" << result << ")" << std::endl;
	}
	
	if (boost::math::isinf(result)) std::cout << "inf: gaussian rand result: " << "(" << result << ")" << std::endl;
		 
	return result;
}

double Utilities::NormalVariate()
{
	//return gen_normal_3(generator);
	return gaussianRand();
}
