//***************************************************************************
//   thermostats.cpp
//   namespace thermostats implementation file
//
//***************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <ctime>
#include <stdexcept>

#include "thermostats.hpp"
#include "rv_library.hpp"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//    Andersen():  Implements Andersen's stochastic thermostat
//
//
//	    Based upon algorithm from H. C. Andersen 1980
//		with modified for use with vortex simulation
// 	 	Jensen 1993.		
//
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void thermostats::Andersen(double T_, double kB_, double eta_, double dt_, double tau_, double * tempForce_)
{
    double p=dt_/tau_;
    double tempforce=0;
	if (p>rv::MT_rand_U())
	{	
		
		tempforce =  sqrt(2*T_*kB_*eta_/dt_/p)*rv::MT_rand_N();
	}
    
	double theta = rv::MT_rand_U()*2*M_PI;
			
	tempForce_[0]=tempforce*sin(theta);
	tempForce_[1]=tempforce*cos(theta);
	
	/*if (tempForce_[0]!=tempForce_[0] || tempForce_[1] != tempForce_[1]) {
		std::cout << "temp nan" << "(" << tempForce_[0] << ", " << tempForce_[1] << ")" << std::endl;
		tempForce_[0]=0;
		tempForce_[1]=0;
	}
	
	if (boost::math::isinf(tempForce_[0]) || boost::math::isinf(tempForce_[1]))
	std::cout << "t: " << sim->get_t() << "temperature inf" << "(" << tempForce_[0] << ", " << tempForce_[1] << ")" << std::endl;
	*/
}

//***************************************************************************
//
//  End
//
//***************************************************************************

 
