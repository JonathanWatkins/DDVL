//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  GeometryShearedWall.h
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#ifndef GeometryShearedWallHPP
#define GeometryShearedWallHPP
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "GeometryBase.hpp"

#include <list>
#include <string>

#include "CParticle.hpp"
#include "CDelLine.hpp"


class CSimulation;
class FileOutput;
class BinnedAccumulator;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  class GeometryShearedWall
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class GeometryShearedWall : public GeometryBase
{
    public:
        GeometryShearedWall(CSimulation * sim_);
        ~GeometryShearedWall();
        
        
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
		
		void CheckEscapedVortices();
        void InitialiseVortices();
		void WrapVorticesX(std::list<CParticle>& jList);
		void WrapVorticesY(std::list<CParticle>& jList);
		void WrapVorticesXY(std::list<CParticle>& jList);
		void DoWrapX(std::list<CParticle>::iterator p, std::list<CParticle>& jList);
		void DoWrapY(std::list<CParticle>::iterator p, std::list<CParticle>& jList);
		void DoWrapXY(std::list<CParticle>::iterator p, std::list<CParticle>& jList);
		void TestX(std::list<CParticle>::iterator p);
		void TestY(std::list<CParticle>::iterator p);
		void KeepParticlesInSimBoxX();
		void KeepParticlesInSimBoxY();
		void KeepParticlesInSimBoxXY();
		void UserUpdates();
		
		void InitialiseFiles();
		void LoadBatchFile();
		void InitialiseParameters();
		void GetPeriodicity();
				 
    	// Analysis functions
		void OutputFinalParticlePositions();
		void OutputParticlePositions();
		void OutputAverages();
		void OutputParticleCount();
		void OutputVxofyProfile();
		void CalculateVxofyProfile();
		void OutputVxofyEvolveProfile();
		
		
		// User Defined functions
		void MoveTopCE();

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// class specific variables
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX		
        
        CSimulation * sim;
        
        std::list<CParticle> * AParticlesList;
        std::list<CParticle> * OtherParticlesList;
        std::list<CParticle> * triangulatedParticlesList;
        std::list<CDelLine> * triangulatedLinesList;
        
		double Phi;
		double a0;
		double b0;
		double dt;
		double forcerange;
		std::string pos_file_name;
        
        double binsize;
        
        bool wrapx;
        bool wrapy;
        double delx;
		double dely;
	
        
        double xlo,ylo, xhi,yhi;
        
        FileOutput * fout;
		
		BinnedAccumulator * Vxofy;
		
        // temp variables for sheared wall jobs
        double topwallvel;
		
		std::string integratorType;
    
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
