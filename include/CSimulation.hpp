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
#include "GeometryBase.hpp"

class CParallelEulerIntegrator;


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
	
	int get_geometry() const {	return geometry; }
	
	std::ofstream* get_FS(const std::string & filestr) {return fileOutputter.getFS(filestr);}
	
	int get_t() { return t; }  // sim owns

	int Geta0() {return geom->Geta0(); }

	int get_simulation_time() { return simulation_time; }  // sim owns
		
	std::string get_jobBatchFileLocation() const { return jobBatchFileLocation; }  // sim owns
	
	double get_xlo()  { return geom->GetXLo(); }  

	double get_xhi() { return geom->GetXHi(); }

	double get_ylo()  { return geom->GetYLo(); }
	
	double get_yhi() { return geom->GetYHi(); }
	
	double get_M2Average() const;  // integrator owns
	
	double get_time() const; // integrator owns
	
	double get_M2FullAverage() const; // integrator owns
	
	GeometryBase * get_geom() { return geom; }
	
	int get_triangulationInterval() {return triangulationInterval; }
	int get_framedataInterval() {return framedataInterval; }
	
	template <class classA>
	void ReadVariableFromBatchFile(classA & to, const std::string & var_);	
	
	
	
private:
	
	void InitialiseFiles();
	
	void DoStep();
	
	void DelaunayTriangulation();
	
	void CalculateAndOutputFinishTime();
	
	void CopyJobBatchFile();
	
	void ReadGeometryType();
	
	void ConfigureSimulation();
	
	void OutputSimulationTimes();
	
	void AssignJobNumber();
	
	GeometryBase * CreateGeometry();
	
	GeometryBase * geom;
	
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
	int framedataInterval;
	int triangulationInterval;
	
	// timing variables
	int DTcount;
	int fcount;
	double DTtime;
	double ftime;
	
	
	int geometry;
	
	CFileOutput fileOutputter;
	
};

template <class classA>
void CSimulation::ReadVariableFromBatchFile(classA & to, const std::string & var_)
{
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(jobBatchFileLocation, pt);
	
	// geometry variables
	to =  pt.get<classA>(var_);
	
}




#endif

