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


// 2D
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
	
}

// 3D
void thermostats::Andersen3D(double T_, double kB_, double eta_, double dt_, double tau_, double * tempForce_)
{
    double p=dt_/tau_;
    double tempforce=0;
	if (p>rv::MT_rand_U())
	{	
		
		tempforce =  sqrt(2*T_*kB_*eta_/dt_/p)*rv::MT_rand_N();
	}
    
	double theta = rv::MT_rand_U()*M_PI;
	double phi = rv::MT_rand_U()*2*M_PI;
			
	tempForce_[0]=tempforce*sin(theta)*cos(phi);
	tempForce_[1]=tempforce*sin(theta)*sin(phi);
	tempForce_[2]=tempforce*cos(theta);
	
	
}


//***************************************************************************
//
//  End
//
//***************************************************************************

 
