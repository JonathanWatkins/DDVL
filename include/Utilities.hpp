#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

namespace Utilities
{

	/*template<class T>
	double gen_normal_3(T &generator)
	{
		return = generator();
	}
*/

	// Using seed for mt as time(0) does not allow runs to be repeated with same seeds
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> >	
	        mygenerator (boost::mt19937(time(0)), boost::normal_distribution<>());
	
	extern const double pi;

	double gaussianRand();
	
	double NormalVariate();

}

#endif
