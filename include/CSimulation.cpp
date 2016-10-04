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
#include "DTwrapper.hpp"
#include "rv_library.hpp"
#include "FileOutput.hpp"

// IntegratorBase Types
#include "CParallelEulerIntegrator.hpp"
#include "CParallelEulerFMAIntegrator.hpp"
#include "IntegratorBuckledSubstrate.hpp"


// GeometryBase Types
#include "GeometryChannel.hpp"
#include "GeometryTube.hpp"
#include "GeometryCustom.hpp" 
#include "GeometryWedge.hpp" 
#include "GeometryOscWall.hpp" 
#include "GeometryShearedWall.hpp" 
#include "GeometryBuckledSubstrate.hpp"


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
#include <chrono>
#include <ratio>

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
	//MonitorPeriod=0;
	startTime=0;
	endTime=0;
	seedtime = 0;
	
	jobnum="";
	initialised=false;
	geometry="";
	jobtag="";
	t=0;
	framedataInterval=0;
	triangulationInterval=0;
	
	
	DTcount=0;
	fcount=0;
	DTtime=0;
	ftime=0;
	
	fout = new FileOutput;
		
	//CVersion version;



}

CSimulation::~CSimulation()
{

	delete geom;
	delete integrator;
	delete fout;
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
	if (0==t%framedataInterval) std::cout << "t: " << t << std::endl;
	
	CalculateAndOutputFinishTime();
	geom->PerStepUpdates();
	
	clock_t startclock = clock();
	integrator->Integrate();
	ftime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	fcount++;
	
	startclock = clock();
	DelaunayTriangulation();
	//delVortexList=vorticesList;	
	DTtime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	DTcount++;
	
}

int CSimulation::Initialise(std::string jobBatchFileLocation_)
{
	// set initial variables to start values for the simulation
	//version.set_versionStr("1.0.0");	
	startTime=clock();
	seedtime = time(0);
	lasttime=std::chrono::steady_clock::now();
	jobBatchFileLocation=jobBatchFileLocation_;
		
	srand ( seedtime );
		
	
	ReadVariableFromBatchFile<std::string>(geometry,"Header.geometry");  // reads geometry from batch file
		
	ReadVariableFromBatchFile(simulation_time,"Header.simulationTime");
	
	ReadVariableFromBatchFile<std::string>(jobtag, "Job.jobtag");  
	
	ReadVariableFromBatchFile(triangulationInterval, "GeneralParameters.triangulationInterval");
	
	ReadVariableFromBatchFile(framedataInterval, "GeneralParameters.framedataInterval");
	
	AssignJobNumber();
	
	InitialiseFileOutput(); // must be done before initialise geometry
		
	geom = CreateGeometry();  
	
	geom->InitialiseGeometry();
		
	integrator = SelectIntegrator(); // geometry knows the integrator to use
	
	integrator->Initialise();
	
	CopyJobBatchFile();
	
	initialised=true;
	
	std::cout << "Simulation initialised.\n\n";
		
	return 0;
}

void CSimulation::InitialiseFileOutput()
{
	
	// make new directory for data
	std::cout << "Initialising FileOutput..." << std::endl;
	fout->SetJobDirectory(jobnum);
	
	
	
	std::cout << "   FileOutput Initialised.\n\n";
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
	std::cout << "All data written to " << jobnum << std::endl;
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
		std::chrono::steady_clock::time_point newtime = std::chrono::steady_clock::now();
		//MonitorPeriod=newtime-lasttime;
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(newtime - lasttime);
		lasttime=newtime;
		
		double MonitorPeriod = time_span.count();

		double hours = floor(MonitorPeriod*(simulation_time-t)/1000.0/60.0/60.0);
		double minutes = (int)(MonitorPeriod*(simulation_time-t)/1000.0/60.0)%60;
		double seconds = (int)(MonitorPeriod*(simulation_time-t)/1000.0)%60;
		std::cout << "Estimated finish time : " << hours << "h " << minutes << "m " << seconds << "s" <<  std::endl;
		
	}
	
}


GeometryBase * CSimulation::CreateGeometry()
{
    //char o_type = inp_->GetOtype();
    
    if(geometry=="channel") return new GeometryChannel(this);
    if(geometry=="tube") return new GeometryTube(this);
    if(geometry=="custom") return new GeometryCustom(this);
    if(geometry=="wedge") return new GeometryWedge(this);
    if(geometry=="oscwall") return new GeometryOscWall(this);
	if(geometry=="shearedwall") return new GeometryShearedWall(this);
	if(geometry=="buckledsubstrate") return new GeometryBuckledSubstrate(this);
	
	
	std::stringstream oss;
	oss << "CSimulation:CreateGeometry()  Bad geometry selected. " << geometry;  
	throw std::runtime_error(oss.str());
	
}

IntegratorBase * CSimulation::SelectIntegrator()
{
	
	std::string integrator_type = geom->GetIntegratorType(); 
	
    if(integrator_type=="ParallelEulerIntegrator") return new CParallelEulerIntegrator(this);
    if(integrator_type=="ParallelEulerFMAIntegrator") return new CParallelEulerFMAIntegrator(this);
    if(integrator_type=="IntegratorBuckledSubstrate") return new IntegratorBuckledSubstrate(this);
    
    

	std::stringstream oss;
	oss << "CSimulation:SelectIntegrator()  Bad integrator selected. " << integrator_type;  
	throw std::runtime_error(oss.str());
	
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
	
	if (t%triangulationInterval!=0) return;
	std::list<CParticle> tmp;
	geom->AddParticlesForDT(tmp);
	
	std::list<CParticle> * trilist = geom->GetTriangulatedParticlesList();
	trilist->clear();
	
	std::list<CParticle>::iterator it = trilist->end();
	
	trilist->insert(it,tmp.begin(),tmp.end());
	
	//std::cout << "num particles: " << tmp.size() << std::endl;
	//return;
	
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	// Temp code for testing the wrapping in the tube geomery
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	/*if (t==285)
	{
		fout->AddFileStream("wraptest","wraptest.txt");
		std::stringstream oss;
		for (std::list<CParticle>::iterator p = tmp.begin();
			p != tmp.end(); ++p)
		{
			oss << p->get_type() << " " << p->get_x() << " " << p->get_y() << std::endl;
			
		}
		
		
		fout->RegisterOutput("wraptest", oss.str());
	}*/
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	// End of Temp code for testing the wrapping in the tube geomery
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
	
	//ComputationalGeometry::DelaunayTriangulation(tmp,geom->GetTriangulatedParticlesList(),geom->GetTriangulatedLinesList());	
}

double CSimulation::get_M2Average() { return integrator->GetM2Average(); }
	
double CSimulation::get_time() const { return integrator->Getdt()*t; };  
	
double CSimulation::get_M2FullAverage() const {	return integrator->GetM2FullAverage(); }

double CSimulation::get_forcerange() { return integrator->GetForceRange(); }

double CSimulation::get_dt() { return integrator->Getdt(); }  // sim owns

	
