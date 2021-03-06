//***************************************************************************
//   rv_library.cpp
//   namespace RandomVariableStatisticalFunctions implementation file
//
//***************************************************************************


#include <random>
#include <chrono>

#include "rv_library.hpp"

// Wrapper for legacy support in FMA integrator

double RandomVariableStatisticalFunctions::GetNormalVariate()
{
    return MT_rand_N();
}

// Implements the Mersenne Twister algorithm or U(0,1)

double RandomVariableStatisticalFunctions::MT_rand_U()
{
	static bool NOT_FIRST_TIME = false;
	static std::mt19937 g1;
	if (NOT_FIRST_TIME == false)
	{
		// obtain a seed from the system clock:
  		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		  g1.seed(seed);  // mt19937 is a standard mersenne_twister_engine
  		NOT_FIRST_TIME = true;
	}
	
	
	std::uniform_real_distribution<double> unidist(0.0,1.0);
  
	return unidist(g1); // returns a random number U(0,1)
	
	
}

double RandomVariableStatisticalFunctions::MT_rand_N()
{
	static bool NOT_FIRST_TIME = false;
	static std::mt19937 g1;
	if (NOT_FIRST_TIME == false)
	{
		// obtain a seed from the system clock:
  		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		  g1.seed(seed);  // mt19937 is a standard mersenne_twister_engine
  		NOT_FIRST_TIME = true;
	}
	
	
	std::normal_distribution<double> normdist(0.0,1.0);
  
	return normdist(g1); // returns a random number N(0,1)
	
	
}

uint_fast32_t RandomVariableStatisticalFunctions::MT_rand()
{
	static bool NOT_FIRST_TIME = false;
	static std::mt19937 g1;
	if (NOT_FIRST_TIME == false)
	{
		// obtain a seed from the system clock:
  		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		  g1.seed(seed);  // mt19937 is a standard mersenne_twister_engine
  		NOT_FIRST_TIME = true;
	}
	
	return g1(); // returns a random number
	
	
}


//***************************************************************************
//
//  End
//
//***************************************************************************

 
