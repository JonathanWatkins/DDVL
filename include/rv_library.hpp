//***************************************************************************
//  rv_library.hpp
//***************************************************************************
#ifndef rv_libraryH
#define rv_libraryH
//***************************************************************************
#include <stdint.h>

namespace RandomVariableStatisticalFunctions
{
   

    //***************************************************************************
    //  Prototypes
    //***************************************************************************
		
		uint_fast32_t MT_rand();						//Random positive integers (32 bit) using Mersenne Twister and c++11
		double MT_rand_U();					//U(0,1) distribution using Mersenne Twister and c++11
		double MT_rand_N();					//N(0,1) distribution using Mersenne Twister and c++11
		double GetNormalVariate();  //legacy support wraps MT_rand_N()
}

namespace rv = RandomVariableStatisticalFunctions;    //alias


#endif 
//***************************************************************************
//  End
//***************************************************************************
