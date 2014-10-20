//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryBase.hpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryBaseHPP
#define GeometryBaseHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <list>
#include <string>

class CParticle;
class CDelLine;


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryBase
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryBase
{
    public:
        virtual ~GeometryBase(){}
        
        virtual void AddParticlesForDT(std::list<CParticle> & vorticesList_) const = 0;
        virtual void PerStepAnalysis() = 0;
		virtual void EndofSimAnalysis() = 0;
		virtual void PerStepUpdates() = 0;
		virtual void InitialiseGeometry() = 0;
		
		virtual void GetIParticles(std::list<CParticle>& vorticesList_) const = 0;  
		virtual void GetJParticles(std::list<CParticle>& vorticesList_) const = 0;
				
		virtual double GetXLo() const = 0;
		virtual double GetXHi() const = 0;
		virtual double GetYLo() const = 0;
		virtual double GetYHi() const = 0;
		virtual double Geta0() const = 0;
		
		
		virtual std::list<CParticle> * GetTriangulatedParticlesList() = 0;
		virtual std::list<CDelLine> * GetTriangulatedLinesList() = 0;
		
		
             
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
