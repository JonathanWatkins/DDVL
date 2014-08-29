//***************************************************************************
//  sca_library.h
//***************************************************************************
#ifndef SCA_LIBRARY_HPP
#define SCA_LIBRARY_HPP
//***************************************************************************

#include <list>

class CParticle;

namespace STLcontainerAlgorithms
{
	
	void ParticlesListCopy(std::list<CParticle> a_, std::list<CParticle> b_, bool copyTrajectories_);
    
}

namespace sca = STLcontainerAlgorithms;    //alias

//***************************************************************************
#endif 
//***************************************************************************
//  End
//***************************************************************************
