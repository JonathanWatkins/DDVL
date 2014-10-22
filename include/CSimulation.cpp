#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

// Class header
#include "CSimulation.hpp"

// Custom classes
#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CDelTriangle.hpp"
#include "CBin.hpp"
#include "CLineIDs.hpp"
#include "CVersion.hpp"
#include "CCell.hpp"
#include "CRunningStats.hpp"
#include "delaunay.hpp"
#include "rv_library.hpp"
#include "CParallelEulerIntegrator.hpp"

// GeometryBase Types
// #include "GeometryChannel.hpp"
#include "GeometryTube.hpp"

//#include "CParameter.hpp"

// STL classes
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <set>
#include <vector>
#include <omp.h>
#include <limits.h>
#include <iomanip>

// Boost libraries
//#include <boost/ptr_container/ptr_list.hpp>
//#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

bool xSort (CParticle lhs_,CParticle rhs_)
{
	return (lhs_.get_x()<rhs_.get_x());
}



CSimulation::CSimulation()
{
	running = false;
	paused = false;
	simulation_time=0;
	MonitorPeriod=0;
	startTime=0;
	endTime=0;
	seedtime = 0;
	lasttime=0;
	std::string jobnum="";
	initialised=false;
	geometry=0;
	jobtag="";
	t=0;
	framedataInterval=0;
	triangulationInterval=0;
	
	
	DTcount=0;
	fcount=0;
	DTtime=0;
	ftime=0;
	
	// encapsulate
	/*M2=0;
	M2Full=0;
	M2Sum=0;
	M2FullSum=0;
	avXVel=0;
	avYVel=0;
	framedataInterval=5;
	cellSize=0;
	binsize=0;
	Nv=0;
	Nmis=0;
	Nd=0;
	thermostat="";
	lorentzForce=0;
	jobtag="";

	Ap=1;
	DTcount=0;
	fcount=0;
	DTtime=0;
	ftime=0;

	applyBathVelocities=false;
	applyStiffBath=false;
	applyBounceBack=false;
	applyMaxVelocities=false;
	alt_pos_file=false;
	pos_file_name="";
	alt_pins_file=false;
	pins_file_name="";
	f0_rcut_correction=0;
	f0bath_rcut_correction=0;
	
	Av=0;
	Rv=0;
	epsilon=0;
	sigma=0;
	vvForce=0;

	*/
	
	

}

void CSimulation::Run()
{
	if (initialised==false)
	{
		throw std::runtime_error("Simulation not initialised.");
	}
	
	std::cout << "Simulation running ..." << std::endl; 
	std::cout << "-------------------------------------------------" << std::endl << std::endl;
	
	for (t=1; t<=simulation_time; ++t)
	{
		DoStep();
		geom->PerStepAnalysis();

	}
	std::cout << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
	std::cout << "Simulation finished." << std::endl << std::endl; 
	geom->EndofSimAnalysis();
	OutputSimulationTimes();
}

void CSimulation::DoStep()
{
	if (0==t%100) std::cout << "t: " << t << std::endl;
	CalculateAndOutputFinishTime();
	geom->PerStepUpdates();
	
	clock_t startclock = clock();
	integrator->Integrate();
	ftime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	fcount++;
	
	startclock = clock();
	//DelaunayTriangulation();
	//delVortexList=vorticesList;	
	DTtime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	DTcount++;
	
}

int CSimulation::Initialise(std::string jobBatchFileLocation_)
{
	// set initial variables to start values for the simulation
	version.set_versionStr("1.0.0");	
	startTime=clock();
	seedtime = time(0);
	lasttime=time(0);
	jobBatchFileLocation=jobBatchFileLocation_;
		
	srand ( seedtime );
		
	// which order should these go XXXXXXXXXXXXXXXXXXXXXXXXXXx
	
	ReadVariableFromBatchFile(geometry,"Header.geometry");  // reads geometry from batch file
	
	geom = CreateGeometry();
	
	geom->InitialiseGeometry();
	
	ReadVariableFromBatchFile(simulation_time,"Header.simulationTime");
	
	ReadVariableFromBatchFile<std::string>(jobtag, "Job.jobtag");  
	
	ReadVariableFromBatchFile(triangulationInterval, "GeneralParameters.triangulationInterval");
	
	ReadVariableFromBatchFile(framedataInterval, "GeneralParameters.framedataInterval");
	
	AssignJobNumber();
		
	// run initial functions
	
	InitialiseFiles();
	
	// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
	// calculate parameters
	
	//A=2*kB*temp/eta;
	
	// initialise integrator
	
	integrator = new CParallelEulerIntegrator(this);
	
	integrator->Initialise();
	
	CopyJobBatchFile();
	
	initialised=true;
	
	std::cout << "Simulation initialised.\n\n";
		
	return 0;
}

void CSimulation::InitialiseFiles()
{
	
	// make new directory for data
	std::cout << "Initialising files..." << std::endl;
	fileOutputter.setJobDirectory(jobnum);
	
	// add files to outputter
		
	fileOutputter.addFileStream("jobheader", "jobheader.ini");
	fileOutputter.addFileStream("posfile", "posdata.txt");
	fileOutputter.addFileStream("guifile", "guidata.dat");
	fileOutputter.addFileStream("framevel", "framevel.txt");
	fileOutputter.addFileStream("Nd", "Nd.txt");
	fileOutputter.addFileStream("avfile", "averagesdata.txt");
	fileOutputter.addFileStream("pinsfile", "pinsdata.txt");
	fileOutputter.addFileStream("wraptest", "wraptest.txt");
    
	std::cout << "   Files initialised.\n\n";
}

void CSimulation::OutputSimulationTimes()
{
	std::cout << std::endl << "Simulation times\n----------------\n";
	endTime = clock();
	
	std::cout << "Run Time: " << (endTime-startTime)/(double)CLOCKS_PER_SEC << std::endl;
	
	std::cout << "Total DTtime: " << DTtime << std::endl;
	std::cout << "Total ftime: " << ftime << std::endl;
	
	
	std::cout << "DTtime: " << DTtime/(double)DTcount << " per iteration (" << DTcount << ")" << std::endl;
	std::cout << "ftime: " << ftime/(double)fcount << " per iteration (" << fcount << ")"<< std::endl;
	std::cout << "----------------" << std::endl;
}

void CSimulation::CopyJobBatchFile()
{
	// copy jobheader to job directory and change name to jobheader.ini
	std::ifstream f1(jobBatchFileLocation.c_str(), std::fstream::binary);
  
  // this is a hack
  std::ostringstream oss;
	oss << jobnum << "//jobheader.ini";
  std::ofstream f2(oss.str().c_str(), std::fstream::trunc|std::fstream::binary);
  f2 << f1.rdbuf();
	// this is a hack
 
}

void CSimulation::CalculateAndOutputFinishTime()
{
	
	if (0==t%1000)
	{
		double newtime=time(0);
		MonitorPeriod=newtime-lasttime;
		lasttime=newtime;

		double hours = floor(MonitorPeriod*(simulation_time-t)/1000.0/60.0/60.0);
		double minutes = (int)(MonitorPeriod*(simulation_time-t)/1000.0/60.0)%60;
		double seconds = (int)(MonitorPeriod*(simulation_time-t)/1000.0)%60;
		std::cout << "Estimated finish time : " << hours << "h " << minutes << "m " << seconds << "s" <<  std::endl;
		
	}
	
}


GeometryBase * CSimulation::CreateGeometry()
{
    //char o_type = inp_->GetOtype();
    
    switch(geometry)
    {
        //case 0:  return new GeometryChannel(*this);            break;
        case 1:  return new GeometryTube(this);             break;
        default:   throw std::runtime_error("CSimulation:CreateOption()  Bad character");
    }
}



void CSimulation::AssignJobNumber()
{
	// make jobnum from seedtime and jobtag from jobbatch file	
	std::stringstream oss;
	oss << "job" << seedtime << "-" << jobtag; 
	jobnum=oss.str();
	std::cout << "JobNum: " << jobnum << std::endl << std::endl;
	
}

void CSimulation::DelaunayTriangulation()
{
	
	std::list<CParticle> tmp;
	geom->AddParticlesForDT(tmp);
	ComputationalGeometry::DelaunayTriangulation(tmp,geom->GetTriangulatedParticlesList(),geom->GetTriangulatedLinesList());	
}

double CSimulation::get_M2Average() const { return integrator->GetM2Average(); }
	
double CSimulation::get_time() const { return integrator->Getdt()*t; };  
	
double CSimulation::get_M2FullAverage() const {	return integrator->GetM2FullAverage(); }

double CSimulation::get_forcerange() { return integrator->GetForceRange(); }
	
