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

class GeometryBase;

#define channel		 	0
#define tube 				1

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
	
	double get_cellSize() const { return cellSize; };
	
	double get_forceRange() const { return forceRange; };
	
	double get_time() const { return dt*t; };
	
	double get_epsilon()  { return epsilon;  };
	
	double get_sigma()  { return sigma;  };
	
	double get_frame_force_d() const { return frame_force_d; };	
	
	double get_frame_force_t() const { return frame_force_t; };	
	
	double get_av_force_d() const { return av_force_d.get_mean(); };	
	
	double get_av_force_t() const { return av_force_t.get_mean(); };	
	
	double get_a0() const { 	return a0; }

	double get_b0() const {	return b0; }

	double get_lambda() const { return lambda; }
		
	double get_f0bath() const { return f0bath; }
		
	double get_f0() const { return f0; }
	
	double get_binsize() const { return binsize; }

	double get_disorderRange() const {	return disorderRange; }
		
	double get_disorderStrength() const { return disorderStrength; }

	double get_f0_rcut_correction() const { return f0_rcut_correction; }

	double get_f0bath_rcut_correction() const { return f0bath_rcut_correction; }

	int get_geometry() const {	return geometry; }

	double get_Av() const { return Av; }

	double get_Rv() const { return Rv; }

	std::list<CParticle>* get_vorticesList() {	return &vorticesList; }

	std::list<CParticle>* get_delVortexList() { return &delVortexList; }

	std::list<CParticle>* get_pinsList() {	return &pinsList; }

	std::list<CDelLine>* get_delLinesList() { return &delLinesList; }

	std::list<CParticle>* get_disorderList() {	return &disorderList; }

	int get_t() { return t; }

	int get_simulation_time() { return simulation_time; }
	
	std::string get_jobBatchFileLocation() const { return jobBatchFileLocation; }
	
	int get_Nv() const { return Nv; }

	double get_M2() const { return M2; }

	double get_M2Full() const { return M2Full; }

	double get_temp() const { return temp; }

	double get_A() const { return A; }

	double get_Ap() const { return Ap; }

	double get_dt() const { return dt; }
	
	double get_eta() const { return eta; }
	
	bool get_applyStiffBath() const { return applyStiffBath; }
	
	bool get_applyBathVelocities() const { return applyBathVelocities; }
	
	bool get_applyMaxVelocities() const { return applyMaxVelocities; }
	
	double get_tau() const { return tau; }
	
	double get_kB() const { return kB; }
	
	int get_vvForce() const { return vvForce; }
	
	double get_lorentzForce() const { return lorentzForce; }
	
	double get_Phi() const { return Phi; }
	
	std::string get_thermostat() const { return thermostat; }
	
	std::string GetPosFileName() const {return pos_file_name; }
	
	std::string GetPinsFileName() const {return pins_file_name; }
	
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
	
	void DelaunayTriangulation(std::list<CParticle> vorticesList_);
	
	void CalculateAvVel();
	
	void CalculateFinishTime();
	
	void CopyJobBatchFile();
	
	void ReadJobBatchFile();
	
	void ConfigureSimulation();
	
	void OutputFinalVortexPositions();
	
	void OutputPinsList(); 											// Outputs the list of pins.
	
	void OutputVortexPositions(); 							// Outputs positions and coord num of vortices using delVortexList
	
	void OutputSimulationTimes();
	
	GeometryBase * CreateGeometry();
	
	GeometryBase * geom;
	
	// This is the location of the batch file containing
	// specific run options
	std::string jobBatchFileLocation;
	
	std::list<CParticle> vorticesList;
	std::list<CParticle> pinsList;
	std::list<CParticle> disorderList;
	std::list<CParticle> delVortexList;
	std::list<CDelLine> delLinesList;
	std::list<CDelTriangle> delTrianglesList;

	CParallelEulerIntegrator* integrator;
		
	bool running;
	bool paused;
	CVersion version;
	int simulation_time;
	int MonitorPeriod;
	int lasttime;
	int seedtime;
	double temp;
	std::string jobnum;
	std::string jobtag;
	double lorentzForce;
	bool initialised;
	int geometry;
	double a0;
	double b0;
	double binsize;
	int numBins;
	double M2;
	double M2Full;
	double M2Sum;
	double M2FullSum;
	int kicks;
	double A;
	int Nv;
	int Nmis;
	double rhov;
	double source_rhov;
	double sink_rhov;
	int ratioOfDefects;
	double pi;
	double Phi;
	double forceRange;
	double eta;
	double kB;
	double mu0;
	double lambda;
  double f0;
  double f0bath;
	double Ap;
	double f0_rcut_correction;
	double f0bath_rcut_correction;
	double Av;
	double Rv;
	double epsilon;
	double sigma;
	int t;
	double dt;
	double tau;
	int triangulationInterval;
	int framedataInterval;
	int outputType;
	bool triangulateReadRun;
	clock_t startTime;
	clock_t endTime;
	std::string thermostat;
	int DTcount;
	int fcount;
	double DTtime;
	double ftime;
	clock_t timer;
	double cellSize;
	int lastchangedSource;
	int lastchangedSink;		
	int vvForce;
	bool applyBathVelocities;
	bool applyStiffBath;
	bool applyBounceBack;
	bool applyMaxVelocities;
	double disorderDensity;
	double disorderStrength;
	double disorderRange;
	bool alt_pos_file;
	bool alt_pins_file;
	std::string pos_file_name;		
	std::string pins_file_name;
	double avXVel;
	double avYVel;
	double frame_force_d;
	double frame_force_t;
	CRunningStats av_force_d;
	CRunningStats av_force_t;
	CFileOutput fileOutputter;
	
};

#endif
