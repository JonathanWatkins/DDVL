//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryBase.hpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryBaseHPP
#define GeometryBaseHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <list>
#include <string>

class CParticle;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryBase
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryBase
{
    public:
        virtual ~GeometryBase(){}
        
        virtual void ReplaceEscapedVortices() const = 0;
        virtual void InitialisePins() = 0;
        virtual void InitialiseVortices() const = 0;
        virtual void AddParticlesForDT(std::list<CParticle> & vorticesList_) const = 0;
        virtual void WrapSystem() const = 0; 
        virtual void InitialiseDisorder() const = 0;
        virtual CParticle GetFirstPin() const = 0;
        virtual double GetRemovalSourceX() const = 0;
        virtual double GetRemovalSinkX() const = 0;
        virtual void UpdateBathDensities() const = 0; // Maintains source and sink densities
				virtual bool AddParticleToBath(std::string location_) const = 0;	//		Adds particle to source or sink
				virtual bool RemoveParticleFromBath(std::string location_) const = 0; 	//		Removes particle from source or sink
				virtual void LoadBatchFile() = 0;
				virtual double GetChannelLength() const = 0;
				virtual double GetChannelWidth() const = 0;
				virtual double GetBathLength() const = 0;
				virtual double GetBathWidth() const = 0;
				virtual double GetSourceBField() const = 0;
				virtual double GetSinkBField() const = 0;
				virtual void WrapVortices(std::list<CParticle>& vorticesList_) const = 0;
	
             
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
