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

#define channel		 	0
#define tube 				1
#define periodic 		2
#define wedge				3
#define BSCCO				4

#define GaussianType		0
#define BesselType			1
#define LJType					2
#define BessLogType			3

#define MAXVFIELDBINS 300
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
	
	// velocityField	
	
	CRunningStats ** vfieldBinVx;//[MAXVFIELDBINS][MAXVFIELDBINS]; // for velocityfield vx
	CRunningStats ** vfieldBinVy;//[MAXVFIELDBINS][MAXVFIELDBINS]; // for velocityfield vy
	
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
	bool drawCoordinateGrid;
	bool calcTrajectories;
	bool showParticleTracker;

	std::string jobtag;
	
	// anealing parameters
	bool aneal;
	double anealTmax;
	int anealCycleTime;
	int anealNumCycles;
	double annealing_T;
	int annealing_endtime;
	
	
	// Size of Lorentz force
	double lorentzForce;
	
	bool initialised;
	// geometry
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
	
	double systemLength;
	double systemWidth;
	
	double channelOffset;

	CParticle firstPin;
	
	CParticle viewpoint;

	double etchsourcex0, etchsourcey0, etchsourcex1, etchsourcey1;
	double etchchannelx0, etchchannely0, etchchannelx1, etchchannely1;
	double etchsinkx0, etchsinky0, etchsinkx1, etchsinky1;
	double etchwedgeangle, etchwedgex1, etchwedgex0, etchwedgey0;
	
	double removesourcex,removesinkx, removesourcey0, removesourcey1; 
	double removechannelx0, removechannelx1, removetopchannely, removebottomchannely;
	double removewedgeangle, removewedgex0, removewedgey0, removewedgex1, removewedgey1, removewedgey2;
			
	double reboundy0, reboundy1; 
	double reboundx0, reboundx1; 
	
	double bouncebacky0, bouncebacky1; 
	double bouncebackx0, bouncebackx1; 
	
	
	double urectx0;
	double urecty0;
	double urectx1;
	double urecty1;

	double bulkx0;
	double bulkx1;
	double bulky0;
	double bulky1;

	double dislocationx0;
	double dislocationx1;
	double dislocationy0;
	double dislocationy1;

	int sourceDensity;
	int sinkDensity;
	int channelDensity;

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
	int drawInterval;  // Must be a product of triangulationInterval
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
	double vfieldBinSize;
	bool flat_channel_ends;
	bool reflected_channel_ends;
	int lastchangedSource;
	int lastchangedSink;		
	int vvForce;
	bool annealing;
	int annealing_time;
	double annealing_factor;
	
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
	std::string alt_pos_file_name;
	std::string alt_pins_file_name;
	
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
	
	~CSimulation()
	{
		//clean up arrays
		for(int i = 0; i < MAXVFIELDBINS; ++i)
		{
			delete [] vfieldBinVx[i];
			delete [] vfieldBinVy[i];
		
		}
		delete [] vfieldBinVx;
		delete [] vfieldBinVy;
		
		
	};
	
	
	// initilise from Job Header
	int Initialise(std::string jobBatchFileLocation_);
	
	void Run();
	
	void DoStep();
	
	void OutputSimulationTimes();
	
	// getters in header file
	
	double get_cellSize() const { return cellSize; };
	
	double get_bouncebackx0() const { return bouncebackx0; };
	
	double get_bouncebacky0() const { return bouncebacky0; };
	
	double get_bouncebackx1() const { return bouncebackx1; };
	
	double get_bouncebacky1() const { return bouncebacky1; };
	
	bool get_applyBounceBack() const { return applyBounceBack; };	
	
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
	
	CParticle get_viewpoint() const {return viewpoint; };
	
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
		
	double get_systemLength() const { return systemLength; }

	double get_systemWidth() const {	return systemWidth; }

	double get_Av() const { return Av; }

	double get_Rv() const { return Rv; }

	std::list<CParticle>* get_vorticesList() {	return &vorticesList; }

	std::list<CParticle>* get_delVortexList() { return &delVortexList; }

	std::list<CParticle>* get_pinsList() {	return &pinsList; }

	std::list<CDelLine>* get_delLinesList() { return &delLinesList; }

	std::list<CParticle>* get_disorderList() {	return &disorderList; }

	int get_t() { return t; }

	int get_simulation_time() { return simulation_time; }

	CParticle get_firstPin() const {	return firstPin; }

	bool is_running() const { return running; }

	bool is_initialised() const { return simulation_initialised; }

	std::string get_finishTimeStr() const { return finishTimeStr; }

	int get_Nv() const { return Nv; }

	std::string get_normaliseSourceStr() const { return normaliseSourceStr; }

	std::string get_normaliseSinkStr() const { return normaliseSinkStr; }

	std::string get_BfieldStr() const { return BfieldStr; }

	bool is_paused() const { return paused; }

	double get_M2() const { return M2; }

	double get_M2Full() const { return M2Full; }

	double get_temp() const { return temp; }

	double get_A() const { return A; }

	double get_Ap() const { return Ap; }

	bool get_drawCoordinateGrid() const { return drawCoordinateGrid; }

	bool get_showParticleTracker() const { return showParticleTracker; }

	int get_runtype() const { return runtype; }

	double get_zoom() const { return zoom; }
	
	double get_dt() const { return dt; }
	
	double get_eta() const { return eta; }
	
	bool get_reflected_channel_ends() const { return reflected_channel_ends; }
	
	bool get_flat_channel_ends() const { return flat_channel_ends; }
	
	bool get_applyStiffBath() const { return applyStiffBath; }
	
	bool get_applyBathVelocities() const { return applyBathVelocities; }
	
	bool get_applyMaxVelocities() const { return applyMaxVelocities; }
	
	double get_tau() const { return tau; }
	
	double get_kB() const { return kB; }
	
	int get_vvForce() const { return vvForce; }
	
	double get_lorentzForce() const { return lorentzForce; }
	
	std::string CSimulation::get_thermostat() const { return thermostat; }
	
	double CSimulation::get_removesourcey0() const { return removesourcey0;	}

	double CSimulation::get_removesourcey1() const { return removesourcey1;	}

	double CSimulation::get_removewedgex0() const { return removewedgex0; }

	double CSimulation::get_removewedgey0() const { return removewedgey0; }

	double CSimulation::get_removewedgex1() const { return removewedgex1; }

	double CSimulation::get_removewedgey1() const { return removewedgey1; }

	double CSimulation::get_removewedgey2() const { return removewedgey2; }

	double CSimulation::get_reboundy0() const { return reboundy0;	}

	double CSimulation::get_reboundy1() const { return reboundy1;	}

	double CSimulation::get_reboundx0() const {	return reboundx0;	}

	double CSimulation::get_reboundx1() const {	return reboundx1;	}
	
				
	//
	//		returns the time average so far of M2
	//    calculated using tempSum[0]*tempSum[0] 
	//
	double get_M2Average() const;
	
	//
	//		returns the time average so far of M2
	//    calculated using (vel_x*dt)^2 
	//
	double get_M2FullAverage() const;
		
	double get_tAvSAvVelX() const;
		
	double get_tAvSAvVelY() const;
	
	//double get_tAvSAvVelXScaled() const;
		
	//double get_tAvSAvVelYScaled() const;
	
	bool get_drawSixFold() const;
	
	bool get_draw() const; // returns true is draw this step. 
	
	void pause_simulation();
	
	void unpause_simulation();
		
	void end_simulation();
	
	void next_t();
	
	void zoom_in();
	
	void zoom_out();
	
	void pan_left();
	
	void pan_right();
	
	void calcBfield();

private:

	//int loop();
	
	void OutputResults();
	
	void calculateForces_old();
	
	void check_for_start();

	void check_for_end();

	void initialise_files();
	
	void initialisePins();
	
	void initialisePinsTube();
	
	void initialisePinsPeriodic();
	
	void initialiseVortices();
	
	void initialiseVorticesTube();
	
	void initialiseVorticesPeriodic();
	
	void initialiseChannelDisorder();
	
	void delaunayTriangulation(std::list<CParticle> vorticesList_);
	
	void removeEscapedVortices();
	
	void removeEscapedVorticesTube(); 
	
	void removeEscapedVorticesPeriodic(); 
	
	void calculateAvVel();
	
	void OutputFinalVortexPositions();
	
	void read_PinsList();
	
	void readSingleDataStep();
	
	void calculateFinishTime();
	
	bool inSystem(double x_, double y_);
	
	// Writes jobBatch information to jobheader 
	// and contains all information required to correctly view the
	// the data using the read run type
	void iniwrite_jobheader();
	
	// reads the run specific ini file
	void iniread_JobBatchFile();
	
		// reads the file written by iniwrite_jobheader()
	void iniread_JobHeader();
	
	void configure_simulation();
	
	// force functions
	
	int what_ivfieldBin(const CParticle &a_) const;
	
	int what_jvfieldBin(const CParticle &a_) const;
	
	double calcSinkB();
		
	double calcSourceB();	
	
	CParticle findClosestParticle (CParticle a_);
	
	CParticle findClosestPin (CParticle a_);
	
	void OutputPinsList(); 											// Outputs the list of pins.
	
	void OutputVortexPositions(); 							// Outputs positions and coord num of vortices using delVortexList
	
	void UpdateTrajectories();									// Update trajectories and calls output at end of simulation
	
	void OutputTrajectory
			(std::list<CParticle>::iterator p_); 		// Outputs a trajectory if vortex escapes (or destroyed)

	int NextTrajectory();												// creates a trajectory and returns the trajectory number
	
	void calculateForces();       							// Calculates forces on particles
	
	
	void ApplyAnnealing();					// Slowly anneals the system

	int FindClosestParticle
		(double x_, double y_, std::vector<CParticle> &test_vector_);			//	Finds the closest particle to a position x, y in a given vector of CParticles, returns index of Particle

	double calcSourceRhoV(); // calculates the denisty of source for LJ and Gaussian simulations
	
	double calcSinkRhoV(); // calculates the denisty of sink for LJ and Gaussian simulations
	
	void SetSimIntervals(); // variable timesteps based on temperature
	
	void updateBathDensities(); // Maintains source and sink densities
	
	bool AddParticleToBath(std::string location_); 	//		Adds particle to source or sink
	
	bool RemoveParticleFromBath(std::string location_); 	//		Removes particle from source or sink
	
	void StartDataCollection();
	
};

#endif
