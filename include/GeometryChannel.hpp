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
#include "CDelLine.hpp"

class CSimulation;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryChannel
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryChannel : public GeometryBase
{
    public:
        GeometryChannel(CSimulation & sim_);
        
        // base functions
				void ReplaceEscapedVortices() const;
        void InitialisePins();
        void InitialiseVortices() const;
        void AddParticlesForDT(std::list<CParticle> & vorticesList_) const;
        void WrapSystem() const;
        void InitialiseDisorder() const;
        CParticle GetFirstPin() const;
        double GetRemovalSourceX() const;
				double GetRemovalSinkX() const;
				void UpdateBathDensities() const;
				bool AddParticleToBath(std::string location_) const;	//		Adds particle to source or sink
				bool RemoveParticleFromBath(std::string location_) const; 	//		Removes particle from source or sink
				void LoadBatchFile();
				void WrapVortices(std::list<CParticle>& vorticesList_) const;
			
				// getters in the header
				double GetChannelLength() const { return channelLength; }
				double GetChannelWidth() const { return channelWidth; }
				double GetBathLength() const { return bathLength; }
				double GetBathWidth() const { return bathWidth; }
				double GetSourceBField() const { return sourceBfield; }
				double GetSinkBField() const { return sinkBfield; }
				
				// class specific functions
				double calcSinkB() const;
				double calcSourceB() const;	
				 
    private:
        
        CSimulation & sim;
        
        std::list<CParticle> * vorticesList;
        std::list<CParticle> * pinsList;
        std::list<CDelLine> * delLinesList;
        
        
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
				double dt;
				double channelOffset;
				std::string pos_file_name;
        std::string pins_file_name;
        std::string jobBatchFileLocation;
        
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
        
        double binsize;
        
 
        
        
              
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
