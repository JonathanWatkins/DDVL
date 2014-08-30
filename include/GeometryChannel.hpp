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

				double etchsourcex0, etchsourcey0, etchsourcex1, etchsourcey1;
				double etchchannelx0, etchchannely0, etchchannelx1, etchchannely1;
				double etchsinkx0, etchsinky0, etchsinkx1, etchsinky1;
        
        int * t;
              
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
