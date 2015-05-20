//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryWedge.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryWedgeHPP
#define GeometryWedgeHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "GeometryBase.hpp"

#include <list>
#include <string>

#include "CParticle.hpp"
#include "CDelLine.hpp"

class CSimulation;
class FileOutput;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryWedge
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryWedge : public GeometryBase
{
    public:
        GeometryWedge(CSimulation * sim_);
        ~GeometryWedge();
        
        
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
		
		
		bool etchParticleWedge(std::list<CParticle>::iterator p_);
		bool removeParticleWedge(std::list<CParticle>::iterator p_);
		
				

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
	
		double etchsourcex0,  
		etchsourcey0,
		etchsourcex1, 
		etchsourcey1;
	
		double etchsinkx0, 
		etchsinky0,
		etchsinkx1, 
		etchsinky1;
	
		double etchchannelx0,
		etchchannely0,
		etchchannelx1,
		etchchannely1;
	
		double removechannelx0,
		removechannelx1,
		removesourcey0,
		removesourcey1;	
		
		double etchwedgex0, etchwedgex1, etchwedgey0, etchwedgeangle;
	
		double removewedgex0, removewedgex1, removewedgey0, removewedgeangle;
			  
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
		              
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
