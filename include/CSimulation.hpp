
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
#include "CCoord.hpp"
#include "CPin.hpp"
#include "CCell.hpp"
#include "CVersion.hpp"
#include "CRowCount.hpp"
#include "CBin.hpp"
#include "CRunningStats.hpp"
#include "CTrajectory.hpp"
#include "CPositionBase.hpp"
#include "CFileOutput.hpp"
#include "CParallelEulerIntegrator.hpp"

class GeometryBase;

#define channel		 	0
#define tube 				1
#define periodic 		2
#define wedge				3
#define BSCCO				4

#define GaussianType		0
#define BesselType			1
#define LJType					2
#define BessLogType			3

#define MAXLINKEDLISTSIZE 100

struct SMoves
{
	double move_x;
	double move_y;
};


class CSimulation
{
	
	friend class CParallelEulerIntegrator;

private:
	
	GeometryBase * geom;
	
	// This is the location of the batch file containing
	// specific run options
	std::string jobBatchFileLocation;
	
	std::list<CParticle> vorticesList;
	std::list<CParticle> pinsList;
	std::list<CParticle> disorderList;
	
	// Triangulation results lists
	std::list<CParticle> delVortexList;
	std::list<CDelLine> delLinesList;
	std::list<CDelTriangle> delTrianglesList;
	
	std::vector<CDelLine> rowCountLinesVector;
	
	std::vector<CRowCount> rowCount;
	
	std::vector<CTrajectory> particleTrajectories;
	std::vector<CParticle> burgers_circuit;
	
	CParallelEulerIntegrator* integrator;
	
	
	// Object status
	bool running;
	bool paused;
	bool simulation_initialised;
	CVersion version;
	
	// initial simulation parameters
	
	int simulation_time;
	int MonitorPeriod;
	int lasttime;
	int seedtime;
	double temp;
	std::string jobnum;
	bool calcTrajectories;
	
	std::string jobtag;
	
	// Size of Lorentz force
	double lorentzForce;
	
	bool initialised;
	int geometry;
	
	double a0;
	double b0;
		
	double sourceBfield;
	double sinkBfield;
	double Bfield;
	double channelLength;
	double channelWidth;
	double bathLength;
	double bathWidth;
	
	// analyse variables
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
	
	// Physical constants
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
	double original_T;
		
	// simulation variables
	int t;
	double dt;
	double tau;
	int triangulationInterval;
	int framedataInterval;
	int outputType;
	int runtype;
	bool triangulateReadRun;
	int starting_time_step;
	clock_t startTime;
	clock_t endTime;
	std::string thermostat;
	int DTcount;
	int fcount;
	double DTtime;
	double ftime;
	clock_t timer;
	double cellSize;
	bool flat_channel_ends;
	bool reflected_channel_ends;
	int lastchangedSource;
	int lastchangedSink;		
	int vvForce;
	
	bool applyBathVelocities;
	bool applyStiffBath;
	bool applyBounceBack;

	bool applyMaxVelocities;
	// channel disorder 
	double disorderDensity; // num pins per a0^2
	double disorderStrength;
	double disorderRange;
	
	bool alt_pos_file;
	bool alt_pins_file;
	
	std::string pos_file_name;		
	std::string pins_file_name;
	
	
	// output strings
	std::string normaliseSourceStr;
	std::string normaliseSinkStr;
	std::string analysisDataStr;
	std::string calculateBinnedBfieldStr1;
	std::string calculateBinnedBfieldStr2;
	std::string finishTimeStr;
	std::string BfieldStr;
	
	// drawing variables
	double vortexSize;
	double zoom;
	
	//**************************** Averages *****************************
	
	//velocities
	double avXVel;
	double avYVel;

	
	// force sizes
	double frame_force_d;
	double frame_force_t;
	CRunningStats av_force_d;
	CRunningStats av_force_t;
	
	
	//*******************************************************************
	
	CFileOutput fileOutputter;

public:
	
	CSimulation();
	
	~CSimulation() {};
	
	
	// initilise from Job Header
	int Initialise(std::string jobBatchFileLocation_);
	
	void Run();
	
	void DoStep();
	
	void OutputSimulationTimes();
	
	// getters in header file
	
	double get_cellSize() const { return cellSize; };
	
	double get_forceRange() const { return forceRange; };
	
	double get_time() const { return dt*t; };
	
	double get_epsilon()  { return epsilon;  };
	
	double get_sigma()  { return sigma;  };
	
	
	// frame forces per particle
	double get_frame_force_d() const { return frame_force_d; };	
	
	double get_frame_force_t() const { return frame_force_t; };	
	
	//average forces per particle
	double get_av_force_d() const { return av_force_d.get_mean(); };	
	
	double get_av_force_t() const { return av_force_t.get_mean(); };	
	
	double get_a0() const { 	return a0; }

	double get_b0() const {	return b0; }

	double get_lambda() const { return lambda; }
		
	double get_f0bath() const { return f0bath; }
		
	double get_f0() const { return f0; }

	double get_disorderRange() const {	return disorderRange; }
		
	double get_disorderStrength() const { return disorderStrength; }

	double get_f0_rcut_correction() const { return f0_rcut_correction; }

	double get_f0bath_rcut_correction() const { return f0bath_rcut_correction; }

	int get_geometry() const {	return geometry; }

	double get_channelLength() const {	return channelLength; }

	double get_channelWidth() const { return channelWidth; }

	double get_bathLength() const { return bathLength; }
	
	double get_bathWidth() const { return bathWidth; }

	double get_vortexSize() const { return vortexSize; }
		
	double get_Av() const { return Av; }

	double get_Rv() const { return Rv; }

	std::list<CParticle>* get_vorticesList() {	return &vorticesList; }

	std::list<CParticle>* get_delVortexList() { return &delVortexList; }

	std::list<CParticle>* get_pinsList() {	return &pinsList; }

	std::list<CDelLine>* get_delLinesList() { return &delLinesList; }

	std::list<CParticle>* get_disorderList() {	return &disorderList; }

	int get_t() { return t; }

	int get_simulation_time() { return simulation_time; }

	std::string get_finishTimeStr() const { return finishTimeStr; }

	int get_Nv() const { return Nv; }

	std::string get_normaliseSourceStr() const { return normaliseSourceStr; }

	std::string get_normaliseSinkStr() const { return normaliseSinkStr; }

	std::string get_BfieldStr() const { return BfieldStr; }

	double get_M2() const { return M2; }

	double get_M2Full() const { return M2Full; }

	double get_temp() const { return temp; }

	double get_A() const { return A; }

	double get_Ap() const { return Ap; }

	int get_runtype() const { return runtype; }

	double get_zoom() const { return zoom; }
	
	double get_dt() const { return dt; }
	
	double get_eta() const { return eta; }
	
	bool get_applyStiffBath() const { return applyStiffBath; }
	
	bool get_applyBathVelocities() const { return applyBathVelocities; }
	
	bool get_applyMaxVelocities() const { return applyMaxVelocities; }
	
	double get_tau() const { return tau; }
	
	double get_kB() const { return kB; }
	
	int get_vvForce() const { return vvForce; }
	
	double get_lorentzForce() const { return lorentzForce; }
	
	double get_sourceBfield() const { return sourceBfield; }
	
	double get_sinkBfield() const { return sinkBfield; }
	
	double get_Phi() const { return Phi; }
	
	
	std::string get_thermostat() const { return thermostat; }
	
	std::string GetPosFileName() const {return pos_file_name; }
	
	std::string GetPinsFileName() const {return pins_file_name; }
		
	double get_M2Average() const;
	
	double get_M2FullAverage() const;
		
	double get_tAvSAvVelX() const;
		
	double get_tAvSAvVelY() const;
	
	void calcBfield();
	
	CParticle get_firstPin() const;

private:

	void OutputResults();
	
	void initialise_files();
	
	void delaunayTriangulation(std::list<CParticle> vorticesList_);
	
	void calculateAvVel();
	
	void OutputFinalVortexPositions();
	
	void calculateFinishTime();
	
	void iniwrite_jobheader();
	
	void iniread_JobBatchFile();
	
	void configure_simulation();
	
	double calcSinkB();
		
	double calcSourceB();	
	
	CParticle findClosestParticle (CParticle a_);
	
	CParticle findClosestPin (CParticle a_);
	
	void OutputPinsList(); 											// Outputs the list of pins.
	
	void OutputVortexPositions(); 							// Outputs positions and coord num of vortices using delVortexList
	
	int FindClosestParticle
		(double x_, double y_, std::vector<CParticle> &test_vector_);			//	Finds the closest particle to a position x, y in a given vector of CParticles, returns index of Particle

	double calcSourceRhoV(); // calculates the denisty of source for LJ and Gaussian simulations
	
	double calcSinkRhoV(); // calculates the denisty of sink for LJ and Gaussian simulations
	
	void updateBathDensities(); // Maintains source and sink densities
	
	bool AddParticleToBath(std::string location_); 	//		Adds particle to source or sink
	
	bool RemoveParticleFromBath(std::string location_); 	//		Removes particle from source or sink
	
	GeometryBase * CreateGeometry();
	
};

#endif
