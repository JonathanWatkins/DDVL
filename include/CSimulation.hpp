#ifndef CSIMULATION_HPP
#define CSIMULATION_HPP

#include <list>
#include <vector>
#include <fstream>
#include <time.h>

#include <boost/function.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CDelTriangle.hpp"
#include "CCell.hpp"
#include "CVersion.hpp"
#include "CBin.hpp"
#include "CRunningStats.hpp"
#include "CFileOutput.hpp"
#include "CParallelEulerIntegrator.hpp"
#include "GeometryBase.hpp"


#define channel		 		0
#define tube 				1
#define custom				2

#define BesselType			1
#define BessLogType			3


struct SMoves
{
	double move_x;
	double move_y;
};


class CSimulation
{
	
	friend class CParallelEulerIntegrator; // should be using pure virtual base instead
	

public:
	
	CSimulation();
	
	~CSimulation() {};
	
	int Initialise(std::string jobBatchFileLocation_);
	
	void Run();
	
	// getters in header file
	
	double get_cellSize() const { return cellSize; };     // integrator owns
	
	double get_forceRange() const { return forceRange; };  // integrator owns
	
	double get_time() const { return dt*t; };  
	
	double get_epsilon()  { return epsilon;  };  // integrator owns
	
	double get_sigma()  { return sigma;  };  // integrator owns
	
	double get_frame_force_d() const { return frame_force_d; };	  // integrator owns
	
	double get_frame_force_t() const { return frame_force_t; };	  // integrator owns
	
	double get_av_force_d() const { return av_force_d.get_mean(); };	// integrator owns
	
	double get_av_force_t() const { return av_force_t.get_mean(); };	// integrator owns
	
	double get_a0() const { 	return a0; }  // geometry owns

	double get_b0() const {	return b0; }  // geometry owns

	double get_lambda() const { return lambda; }   // integrator owns
		
	double get_f0bath() const { return f0bath; }  // integrator owns
		
	double get_f0() const { return f0; }  // integrator owns
	
	double get_binsize() const { return binsize; } 

	double get_f0_rcut_correction() const { return f0_rcut_correction; }  // integrator owns

	double get_f0bath_rcut_correction() const { return f0bath_rcut_correction; }  // integrator owns

	int get_geometry() const {	return geometry; }

	double get_Av() const { return Av; }  

	double get_Rv() const { return Rv; }

	std::list<CParticle>* get_vorticesList() {	return &vorticesList; }

	std::list<CParticle>* get_delVortexList() { return &delVortexList; }

	std::list<CDelLine>* get_delLinesList() { return &delLinesList; }

	int get_t() { return t; }  // sim owns

	int get_simulation_time() { return simulation_time; }  // sim owns
	
	std::string get_jobBatchFileLocation() const { return jobBatchFileLocation; }  // sim owns
	
	double get_temp() const { return temp; }

	double get_A() const { return A; }

	double get_Ap() const { return Ap; }

	double get_dt() const { return dt; }  // integrator owns
	
	double get_eta() const { return eta; }
	
	bool get_applyMaxVelocities() const { return applyMaxVelocities; }
	
	double get_tau() const { return tau; }  // integrator owns
	
	double get_kB() const { return kB; }  // integrator owns
	
	int get_vvForce() const { return vvForce; }  // integrator owns
	
	double get_lorentzForce() const { return lorentzForce; }  // integrator owns
	
	double get_Phi() const { return Phi; }  // integrator owns
	
	std::string get_thermostat() const { return thermostat; }  // integrator owns
	
	double get_xlo()  { return geom->GetXLo(); }  

	double get_xhi() { return geom->GetXHi(); }

	double get_ylo()  { return geom->GetYLo(); }
	
	double get_yhi() { return geom->GetYHi(); }
	
	double get_M2Average() const { return integrator->GetM2Average(); }

	double get_M2FullAverage() const {	return integrator->GetM2FullAverage(); }
	
	GeometryBase * get_geom() { return geom; }
		
	double get_M2Average() const;
	
	double get_M2FullAverage() const;
		
	double get_tAvSAvVelX() const;
		
	double get_tAvSAvVelY() const;
	
	// from geom class
	
	CParticle get_firstPin() const;
	
	double get_channelLength() const;

	double get_channelWidth() const;

	double get_bathLength() const;
	
	double get_bathWidth() const;

private:
	
	void InitialiseFiles();
	
	void DoStep();
	
	void DelaunayTriangulation();
	
	void CalculateAndOutputFinishTime();
	
	void CalculateAndOutputNd();
	
	void CopyJobBatchFile();
	
	void ReadGeometryType();
	
	void ConfigureSimulation();
	
	void OutputSimulationTimes();
	
	void AssignJobNumber();
	
	template <class classA>
	classA ReadVariableFromBatchFile(const std::string & var_);

	
	GeometryBase * CreateGeometry();
	
	GeometryBase * geom;
	
	// This is the location of the batch file containing
	// specific run options
	std::string jobBatchFileLocation;
	
	
	CParallelEulerIntegrator* integrator;
		
	bool running;
	bool paused;
	CVersion version;
	int simulation_time;
	int MonitorPeriod;
	int lasttime;
	int seedtime;
	double startTime;
	double endTime;
	std::string jobnum;
	bool initialised;
	std::string jobtag;
	int t;
	
	// timing variables
	int DTcount;
	int fcount;
	double DTtime;
	double ftime;
	
	
	int geometry;
	
	CFileOutput fileOutputter;
	
};

#endif
