//***************************************************************************
//   sca_library.cpp
//   namespace STLcontainerAlgorithms implementation file
//
//***************************************************************************

#include <iostream>
#include <algorithm>
#include <list>
#include "CParticle.hpp"

#include "sca_library.hpp"

void STLcontainerAlgorithms::ParticlesListCopy(std::list<CParticle> a_, std::list<CParticle> b_, bool copyTrajectories_)
{
	b_.clear();
	for (std::list<CParticle>::iterator p = a_.begin();
			p != a_.end(); ++p)
	{
			CParticle newParticle;
			newParticle.copyAssign_NoTraj(*p);
			b_.push_back(newParticle);
		
	}
	
	
}

//***************************************************************************
//
//  End
//
//***************************************************************************

 
