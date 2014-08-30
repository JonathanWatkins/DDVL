//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryChannel.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryChannelHPP
#define GeometryChannelHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "GeometryBase.hpp"

#include <list>
#include <string>

#include "CParticle.hpp"

class CSimulation;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryChannel
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryChannel : public GeometryBase
{
    public:
        GeometryChannel(CSimulation & sim_);
        
				void ReplaceEscapedVortices() const;
        void InitialisePins();
        void InitialiseVortices() const;
        void AddParticlesForDT(std::list<CParticle> & vorticesList_) const;
        void WrapSystem() const;
        void InitialiseDisorder() const;
        CParticle GetFirstPin() const;
        double GetRemovalSourceX() const;
				double GetRemovalSinkX() const;
				 
    private:
        
        std::list<CParticle> * vorticesList;
        std::list<CParticle> * pinsList;
        
        CParticle firstPin;
        
        double bathLength;
				double bathWidth;
				double channelLength;
				double channelWidth;
				double sourceBfield;
				double sinkBfield;
				double Phi;
				double a0;
				double b0;
				double channelOffset;
				std::string pos_file_name;
        std::string pins_file_name;
        
        int sourceDensity;
				int sinkDensity;
				int channelDensity;
        
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
