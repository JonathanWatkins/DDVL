//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryTube.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryTubeHPP
#define GeometryTubeHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "GeometryBase.hpp"

#include <list>
#include <string>

#include "CParticle.hpp"
#include "CDelLine.hpp"

class CSimulation;
class FileOutput;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryTube
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryTube : public GeometryBase
{
    public:
        GeometryTube(CSimulation * sim_);
        ~GeometryTube();
        
        
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// base functions
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	    void PerStepAnalysis();
		void EndofSimAnalysis();
		void PerStepUpdates();
		void InitialiseGeometry();
		
		void AddParticlesForDT(std::list<CParticle> & vorticesList_);  // returns particles to be triangulated
		std::list<CParticle> * GetIParticles();  // returns particles to be integrated
		void GetJParticles(std::list<CParticle>& vorticesList_);  // returns particles seen by integrated particles
		 
				
		double GetXLo() const { return xlo; }
		double GetXHi() const { return xhi; }
		double GetYLo() const { return ylo; }
		double GetYHi() const { return yhi; }
		double Geta0() const { return a0; }
		
		
		std::list<CParticle> * GetTriangulatedParticlesList() {return triangulatedParticlesList; }
		std::list<CDelLine> * GetTriangulatedLinesList() {return triangulatedLinesList; }
		
		std::string GetIntegratorType() const {return integratorType; }

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// class specific functions
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX		
private:
		
		void ReplaceEscapedVortices();
        void InitialiseVortices();
        double GetRemovalSourceX() const;
		double GetRemovalSinkX() const;
		void UpdateBathDensities();
		bool AddParticleToBath(std::string location_);	//		Adds particle to source or sink
		bool RemoveParticleFromBath(std::string location_); 	//		Removes particle from source or sink
		
		void WrapVorticesY(std::list<CParticle>& jList);
		void DoWrapY(std::list<CParticle>::iterator p, std::list<CParticle>& jList);
		
		void LoadBatchFile();
		double CalcSinkB() const;
		double CalcSourceB() const;	
		void InitialiseRandomMobileParticles();
		void InitialiseCEParticles();
		void InitialiseParameters();
		void InitialiseFiles();
				 
    	// Analysis functions
		void CalculateAndOutputAvVel();
		void OutputFinalParticlePositions();
		void OutputParticlePositions();
		void OutputAverages();
		void CalculateAndOutputNd();
		
		
		
				

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// class specific variables
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX		
        
        CSimulation * sim;
        
        std::list<CParticle> * AParticlesList;
        std::list<CParticle> * OtherParticlesList;
        std::list<CParticle> * triangulatedParticlesList;
        std::list<CDelLine> * triangulatedLinesList;
        
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
		double forcerange;
		std::string pos_file_name;
        std::string pins_file_name;
        
        int sourceDensity;
		int sinkDensity;
		int channelDensity;
        
        double removesourcex;
        double removesinkx;
        
        double removetopchannely;
        double removebottomchannely;
        
		double etchsourcex, etchsinkx;
			  
        double binsize;
        
        double avXVel;
        double avYVel;
        
        double xlo,ylo, xhi,yhi;
        
        int Nd;
		int Nv;
		int Nmis;
		
		FileOutput * fout;
		
		bool wrapx;
		bool wrapy;
        
		std::string integratorType;      
			  
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
