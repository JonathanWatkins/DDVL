//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryChannel.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryChannelHPP
#define GeometryChannelHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "GeometryBase.hpp"

#include <list>

class CParticle;
class CSimulation;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryChannel
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryChannel : public GeometryBase
{
    public:
        GeometryChannel(CSimulation & sim_);
        
				void RemoveEscapedVortices() const;
        void InitialisePins() const;
        void InitialiseVortices() const;
        void AddParticlesForDT() const;
        void WrapSystem() const;
        
    private:
        
        std::list<CParticle> * vorticesList;
        
        double removesourcex;
        double removesourcey0;
        double removesinkx;
        double removesourcey1;
        double removechannelx0;
        double removechannelx1;
        double removetopchannely;
        double removebottomchannely;
        
        int * t;
              
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
