#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "CParticle.hpp"
#include <list>
#include <cmath>


namespace Utilities
{
	
	extern const double pi;

	void wrapVortices(std::list<CParticle>& vorticesList_, double ysize_, double wrapsize_);
	
	void wrapVorticesPeriodic(std::list<CParticle>& vorticesList_, double xsize_, double ysize_, double wrapsize_);
	
	double gaussianRand();

}

#endif
